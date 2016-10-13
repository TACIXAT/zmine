#include "csolv.h"

// TODO: these two are the same as the other ones 
    // figure out how to compile device code in a third file

// shout out to salix alba, you're a wizard mate
// http://stackoverflow.com/a/39862297/1176872
__device__ uint32_t final_elem_idx_1(uint32_t index) {
    float n = index;
    return (uint32_t)ceil((sqrt(8*n+1)-1)/2);
}

__device__ uint32_t final_elem_idx_2(uint32_t index) {
    float x = final_elem_idx_1(index);
    float t = x * (x-1) / 2;
    return (uint32_t)(index - t);
}

__global__ void final_mark_starts(
        uint8_t *hashes, 
        uint32_t *sort_indices,
        uint32_t *off_map, 
        uint32_t r,
        uint32_t hash_count) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < hash_count) {
        uint32_t t_prev_index = (t_index-1) % hash_count; // wrap around at index 0

        uint32_t index = sort_indices[t_index];
        uint32_t prev_index = sort_indices[t_prev_index]; 

        unsigned char* hash = hashes+index*30*sizeof(unsigned char)+r*3;
        unsigned char* prev_hash = hashes+prev_index*30*sizeof(unsigned char)+r*3;

        uint64_t key = ((uint64_t)hash[0]) << 40 | ((uint64_t)hash[1]) << 32 | hash[2] << 24;
        key |= hash[3] << 16 | hash[4] << 8 | hash[5];

        uint64_t prev_key = ((uint64_t)prev_hash[0]) << 40 | ((uint64_t)prev_hash[1]) << 32 | prev_hash[2] << 24;
        prev_key |= prev_hash[3] << 16 | prev_hash[4] << 8 | prev_hash[5];

        if((key ^ prev_key) != 0) {
            off_map[t_index] = 1;
        }
    }
}

__global__ void final_build_start_end_count(
        uint8_t *hashes, 
        uint32_t *sort_indices,
        uint32_t *off_map,
        uint32_t *start, 
        uint32_t *end, 
        uint32_t *count, 
        uint32_t r,
        uint32_t hash_count) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < hash_count) {
        uint32_t t_prev_index = (t_index-1) % hash_count; // wrap around at index 0

        uint32_t index = sort_indices[t_index];
        uint32_t prev_index = sort_indices[t_prev_index]; 

        unsigned char* hash = hashes+index*30*sizeof(unsigned char)+r*3;
        unsigned char* prev_hash = hashes+prev_index*30*sizeof(unsigned char)+r*3;

        uint64_t key = ((uint64_t)hash[0]) << 40 | ((uint64_t)hash[1]) << 32 | hash[2] << 24;
        key |= hash[3] << 16 | hash[4] << 8 | hash[5];

        uint64_t prev_key = ((uint64_t)prev_hash[0]) << 40 | ((uint64_t)prev_hash[1]) << 32 | prev_hash[2] << 24;
        prev_key |= prev_hash[3] << 16 | prev_hash[4] << 8 | prev_hash[5];

        uint32_t lookup_idx = off_map[t_index];

        if((key ^ prev_key) != 0) {
            start[lookup_idx] = t_index;
            end[lookup_idx] = t_prev_index;
        }

        // unsigned long long int* ull_addr = (unsigned long long int*)&count[lookup_idx];
        atomicAdd(&count[lookup_idx], 1);
    }
}

__global__ void final_build_combination_count(
        uint8_t *hashes, 
        uint32_t *sort_indices,
        uint32_t *off_map,
        uint32_t *start, 
        uint32_t *count,
        uint32_t *combination_count_lookup,
        uint32_t r,
        uint32_t size) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        uint32_t lookup_idx = off_map[t_index];
        uint32_t index = sort_indices[t_index];

        if(index == sort_indices[start[lookup_idx]]) {
            uint32_t key_count = count[lookup_idx];
            uint32_t n = key_count - 1;
            combination_count_lookup[lookup_idx] = n * (n + 1) / 2;
        }
    }

}

__global__ void final_map_index_to_prefix(
        uint8_t *hashes, 
        uint32_t *sort_indices,
        uint32_t *off_map,
        uint32_t *comb_count,
        uint32_t *comb_sum,
        uint32_t *comb_prefix,
        uint32_t r,
        uint32_t size) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        uint32_t index = sort_indices[t_index];
        unsigned char* hash = hashes+index*30*sizeof(unsigned char)+r*3;
        
        uint64_t key = ((uint64_t)hash[0]) << 40 | ((uint64_t)hash[1]) << 32 | hash[2] << 24;
        key |= hash[3] << 16 | hash[4] << 8 | hash[5];

        uint32_t lookup_idx = off_map[t_index];

        uint64_t count = comb_count[lookup_idx];
        uint64_t sum = comb_sum[lookup_idx];
        for(int i=(sum-count); i<sum; i++) {
            comb_prefix[i] = lookup_idx;
        }
    }
}


//TODO: remove end
__global__ void final_calculate_pairs(
        uint32_t size,
        uint32_t *comb_prefix_map,
        uint32_t *sort_indices,
        uint32_t *start, 
        uint32_t *end, 
        uint32_t *count,
        uint32_t *combination_count,
        uint32_t *combination_sum,
        uint32_t *ij_buf,
        uint32_t sum_prev_size) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        uint32_t lookup_idx = comb_prefix_map[t_index];

        uint32_t length = count[lookup_idx];
        uint32_t comb_length = combination_count[lookup_idx];

        uint32_t comb_start = combination_sum[lookup_idx] - combination_count[lookup_idx];
        uint32_t local_index = t_index - comb_start;

        uint32_t lm1 = length - 1;

        uint32_t prefix_start = start[lookup_idx];
        uint32_t elem_t_idx_1 = prefix_start + lm1 - final_elem_idx_1(comb_length - local_index);
        uint32_t elem_t_idx_2 = prefix_start + length - final_elem_idx_2(comb_length - local_index);

        uint32_t elem_idx_1 = sort_indices[elem_t_idx_1];
        uint32_t elem_idx_2 = sort_indices[elem_t_idx_2];

        *(ij_buf+2*sum_prev_size+t_index) = elem_idx_1;
        *(ij_buf+2*sum_prev_size+size+t_index) = elem_idx_2;
    }
}

uint32_t cuda_final_calculate_pairs(
        unsigned char **d_hashes, 
        uint32_t *d_sort_indices, 
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes) 
{
    uint32_t hash_count = buf_sizes[r];
    uint32_t lookup_size = sizeof(uint32_t)*hash_count;

    // populate start, end, and count arrays
    uint32_t block_count = hash_count / BLOCK_SIZE;
    if(hash_count % BLOCK_SIZE != 0)
        block_count += 1;

    uint32_t *d_ones;
    cudaMalloc(&d_ones, lookup_size);
    cudaMemset(d_ones, 0, lookup_size); 

    initialize_to_one<<<block_count, BLOCK_SIZE>>>(d_ones, hash_count);

    // offset map table (hash_count)
    // maps t_index -> lookup index
    uint32_t *d_off_map;
    cudaMalloc(&d_off_map, lookup_size);
    cudaMemset(d_off_map, 0, lookup_size); 

    final_mark_starts<<<block_count, BLOCK_SIZE>>>(
        *d_hashes,
        d_sort_indices,
        d_off_map,
        r,
        hash_count);

    thrust::device_vector<uint32_t> dv_off_map(d_off_map, d_off_map+hash_count);
    thrust::device_vector<uint32_t> dv_ones(d_ones, d_ones+hash_count);

    thrust::inclusive_scan(
        dv_off_map.begin(), dv_off_map.end(), dv_off_map.begin());

    thrust::transform(
        dv_off_map.begin(), 
        dv_off_map.end(), 
        dv_ones.begin(), 
        dv_off_map.begin(), 
        thrust::minus<uint32_t>());


    std::vector<uint32_t> h_off_map(dv_off_map.size());
    thrust::copy(dv_off_map.begin(), dv_off_map.end(), h_off_map.begin());
    for(int i=1; i<hash_count; i++)
        if(h_off_map[i] != h_off_map[i-1]+1)
            printf("OFF: %x\n", h_off_map[i]);

    cudaMemcpy(d_off_map, thrust::raw_pointer_cast(dv_off_map.data()), lookup_size, cudaMemcpyDeviceToDevice);
    cudaFree(d_ones);
    uint32_t h_map_length;
    cudaMemcpy(&h_map_length, d_off_map+(hash_count-1), sizeof(uint32_t), cudaMemcpyDeviceToHost);
    h_map_length++;

    printf("MAP SIZE: %x\n", h_map_length);

    // start lookup table (h_map_length elements)
    uint32_t *d_start_lookup;
    cudaMalloc(&d_start_lookup, h_map_length*sizeof(uint32_t));
    cudaMemset(d_start_lookup, 0, h_map_length*sizeof(uint32_t)); 

    // end lookup table (h_map_length elements)
    uint32_t *d_end_lookup;
    cudaMalloc(&d_end_lookup, h_map_length*sizeof(uint32_t));
    cudaMemset(d_end_lookup, 0, h_map_length*sizeof(uint32_t)); 

    // prefix count lookup table
    uint32_t *d_count_lookup;
    cudaMalloc(&d_count_lookup, h_map_length*sizeof(uint32_t));  
    cudaMemset(d_count_lookup, 0, h_map_length*sizeof(uint32_t)); 

    final_build_start_end_count<<<block_count, BLOCK_SIZE>>>(
        *d_hashes, 
        d_sort_indices, 
        d_off_map,
        d_start_lookup, 
        d_end_lookup, 
        d_count_lookup, 
        r, 
        hash_count);

    // combination_count lookup table
    uint32_t *d_combination_count_lookup;
    cudaMalloc(&d_combination_count_lookup, h_map_length*sizeof(uint32_t));  
    cudaMemset(d_combination_count_lookup, 0, h_map_length*sizeof(uint32_t)); 

    // k2
    final_build_combination_count<<<block_count, BLOCK_SIZE>>>(
        *d_hashes, 
        d_sort_indices, 
        d_off_map,
        d_start_lookup, 
        d_count_lookup, 
        d_combination_count_lookup, 
        r, 
        hash_count);

    // sum sizes in combination lookup table
    thrust::device_vector<uint32_t> dv_ccl(d_combination_count_lookup, d_combination_count_lookup+h_map_length);

    thrust::inclusive_scan(
        dv_ccl.begin(), dv_ccl.end(), dv_ccl.begin());

    uint32_t *d_combination_sum = thrust::raw_pointer_cast(dv_ccl.data());
    // std::vector<uint32_t> h_combination_sum(h_map_length);
    // thrust::copy(dv_ccl.begin(), dv_ccl.end(), h_combination_sum.begin());

    // for(int i=0; i<20; i++) {
    //     printf("COMB: %x\n", h_combination_sum[i]);
    // }

    // uint32_t h_combination_sum;
    uint32_t size;
    cudaMemcpy(&size, d_combination_sum+(h_map_length-1), 1*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    buf_sizes[r+1] = size;
    // printf("SIZE: %x\n", size);

    printf("ROUND: %x\n", r);
    printf("SIZE: %x\n", size);

    // map prefixes to thread indices
    uint32_t *d_tidx_prefix_map;
    cudaMalloc(&d_tidx_prefix_map, size*sizeof(uint32_t));

    final_map_index_to_prefix<<<block_count, BLOCK_SIZE>>>(
        *d_hashes, 
        d_sort_indices, 
        d_off_map,
        d_combination_count_lookup, 
        d_combination_sum, 
        d_tidx_prefix_map, 
        r,
        hash_count);

    uint32_t size_sum = 0;
    uint32_t sum_prev_size = 0;
    for(int i=0; i<r+1; i++) {
        size_sum += buf_sizes[1+i];
        if(i != r)
            sum_prev_size += buf_sizes[1+i];
    }

    uint32_t *temp_d_ij_buf;
    cudaMalloc(&temp_d_ij_buf, sizeof(uint32_t)*size_sum*2);
    if(r != 0) {
        cudaMemcpy(temp_d_ij_buf, *d_ij_buf, sum_prev_size*sizeof(uint32_t)*2, cudaMemcpyDeviceToDevice);
        cudaFree(*d_ij_buf);
    }

    cudaMalloc(d_ij_buf, sizeof(uint32_t)*size_sum*2);
    cudaMemcpy(*d_ij_buf, temp_d_ij_buf, sizeof(uint32_t)*size_sum*2, cudaMemcpyDeviceToDevice);
    cudaFree(temp_d_ij_buf);

    block_count = size / BLOCK_SIZE;
    if(size % BLOCK_SIZE != 0)
        block_count += 1;

    // each thread calculates its pair and stores 
    final_calculate_pairs<<<block_count, BLOCK_SIZE>>>(
        size,
        d_tidx_prefix_map,
        d_sort_indices,
        d_start_lookup, 
        d_end_lookup, 
        d_count_lookup,
        d_combination_count_lookup,
        d_combination_sum,
        *d_ij_buf,
        sum_prev_size);

    // free all the crap used to calculate pairs
    cudaFree(d_tidx_prefix_map);
    cudaFree(d_start_lookup);
    cudaFree(d_end_lookup);
    cudaFree(d_count_lookup);
    cudaFree(d_combination_count_lookup);

    return size;
}

void cuda_final_collision_step(
        unsigned char **d_hashes, 
        uint32_t *d_sort_indices, 
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes) 
{
    // check collision on 6 low bytes (uint64_t)
    // create pairs
    uint32_t size = cuda_final_calculate_pairs(d_hashes, d_sort_indices, d_ij_buf, r, buf_sizes);

    uint32_t block_count = size / BLOCK_SIZE;
    if(size % BLOCK_SIZE != 0)
        block_count += 1;

    uint32_t size_sum = 0;
    uint32_t sum_prev_size = 0;
    for(int i=0; i<r+1; i++) {
        size_sum += buf_sizes[1+i];
        if(i != r)
            sum_prev_size += buf_sizes[1+i];
    }
    printf("BLOCK COUNT: %x\n", block_count);
    // check indices
    size = cuda_reduce_pairs(d_ij_buf, r, buf_sizes, size, sum_prev_size, block_count);
    printf("FINAL SIZE: %x\n", size);
}