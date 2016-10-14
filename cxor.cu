// Copyright Douglas Goddard 2016
// Licensed under the MIT license

#include "csolv.h"

// shout out to salix alba, you're a wizard mate
// http://stackoverflow.com/a/39862297/1176872
__device__ uint32_t elem_idx_1(uint32_t index) {
    float n = index;
    return (uint32_t)ceil((sqrt(8*n+1)-1)/2);
}

__device__ uint32_t elem_idx_2(uint32_t index) {
    float x = elem_idx_1(index);
    float t = x * (x-1) / 2;
    return (uint32_t)(index - t);
}

__global__ void build_start_end_count(
        uint8_t *hashes, 
        uint32_t *sort_indices,
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

        uint32_t key = hash[0] << 16 | hash[1] << 8 | hash[2];
        uint32_t prev_key = prev_hash[0] << 16 | prev_hash[1] << 8 | prev_hash[2];

        if((key ^ prev_key) != 0) {
            start[key] = t_index;
            end[prev_key] = t_prev_index;
        }

        atomicAdd(&count[key], 1);
    }
}

__global__ void build_combination_count(
        uint8_t *hashes, 
        uint32_t *sort_indices,
        uint32_t *start, 
        uint32_t *count,
        uint32_t *combination_count_lookup,
        uint32_t r,
        uint32_t size) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        uint32_t index = sort_indices[t_index];
        unsigned char* hash = hashes+index*30*sizeof(unsigned char)+r*3;
        uint32_t key = hash[0] << 16 | hash[1] << 8 | hash[2];

        if(index == sort_indices[start[key]]) {
            uint32_t key_count = count[key];
            uint32_t n = key_count - 1;
            combination_count_lookup[key] = n * (n + 1) / 2;
        }
    }

}

__global__ void map_index_to_prefix(
        uint8_t *hashes, 
        uint32_t *sort_indices,
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
        uint32_t key = hash[0] << 16 | hash[1] << 8 | hash[2];

        uint32_t count = comb_count[key];
        uint32_t sum = comb_sum[key];
        for(int i=(sum-count); i<sum; i++) {
            comb_prefix[i] = key;
        }
    }
}

__global__ void calculate_pairs(
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
        uint32_t prefix = comb_prefix_map[t_index];

        uint32_t length = count[prefix];
        uint32_t comb_length = combination_count[prefix];

        uint32_t comb_start = combination_sum[prefix] - combination_count[prefix];
        uint32_t local_index = t_index - comb_start;

        uint32_t lm1 = length - 1;

        uint32_t prefix_start = start[prefix];
        uint32_t elem_t_idx_1 = prefix_start + lm1 - elem_idx_1(comb_length - local_index);
        uint32_t elem_t_idx_2 = prefix_start + length - elem_idx_2(comb_length - local_index);

        uint32_t elem_idx_1 = sort_indices[elem_t_idx_1];
        uint32_t elem_idx_2 = sort_indices[elem_t_idx_2];

        *(ij_buf+2*sum_prev_size+t_index) = elem_idx_1;
        *(ij_buf+2*sum_prev_size+size+t_index) = elem_idx_2;
    }
}

// TODO: implement merge sort on recovered indices
__global__ void get_reduction(
        //clock_t * clocks,
        uint32_t *reduction, 
        uint32_t r,
        uint32_t *size_buf,
        uint32_t size,
        uint32_t sum_prev_size,
        uint32_t *ij_buf) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        // clock_t start;
        //if(t_index == 0) {
        //     start = clock();
        //}
        uint32_t root_size = (1<<(r+1));
        // uint32_t *roots = (uint32_t *)malloc(sizeof(uint32_t)*root_size);
        uint32_t roots[512];

        // curr = root
        uint32_t curr_node = t_index;
        // depth = round
        uint32_t curr_depth = r;
        // visited = false
        bool visited = false;

        uint32_t curr_size = size;
        uint32_t curr_base = sum_prev_size;

        // node stack
        uint32_t node_stack[8];
        // depth stack
        uint32_t depth_stack[8];
        // size stack
        uint32_t base_stack[8];
        
        int stack_ptr = 0;
        int root_ptr = 0;

        // traverse
        while(root_ptr < root_size) {
            if(!visited) {
                if(!curr_depth) {
                    roots[root_ptr] = curr_node;
                    root_ptr++;
                    visited = true;
                } else {
                    // push
                    node_stack[stack_ptr] = curr_node;
                    depth_stack[stack_ptr] = curr_depth;
                    base_stack[stack_ptr] = curr_base;
                    stack_ptr++;
                    // curr = left
                    curr_node = *(ij_buf + curr_base*2 + curr_node);
                    curr_depth -= 1;
                    curr_size = size_buf[curr_depth];
                    curr_base -= curr_size;
                }
            } else if(visited) {
                curr_node = *(ij_buf + curr_base*2 + curr_size + curr_node);
                if(!curr_depth) {
                    roots[root_ptr] = curr_node;
                    root_ptr++;
                    // pop
                    stack_ptr--;
                    curr_node = node_stack[stack_ptr];
                    curr_depth = depth_stack[stack_ptr];
                    curr_base = base_stack[stack_ptr];
                    curr_size = size_buf[curr_depth];
                    visited = true;
                } else {
                    curr_depth -= 1;
                    curr_size = size_buf[curr_depth];
                    curr_base -= curr_size;
                    visited = false;
                }
            } 
        }

        uint32_t set = 1;
        for(int i=0; i<(root_size/2); i++) {
            for(int j=(root_size/2); j<root_size; j++) {
                if(roots[i] == roots[j]) {
                    set = 0;
                }
            }
        }

        if(set == 0) {
            reduction[t_index] = 0;
        }
    }   
}

__global__ void remap_reduction(
        uint32_t *d_reduction,
        uint32_t *d_mapping,
        uint32_t *old_d_ij_buf,
        uint32_t sum_prev_size,
        uint32_t prev_size,
        uint32_t *new_d_ij_buf,
        uint32_t new_size) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < prev_size) {
        if(d_reduction[t_index]) {
            uint32_t index = d_mapping[t_index];
            uint32_t i = *(old_d_ij_buf+2*sum_prev_size+t_index);
            uint32_t j = *(old_d_ij_buf+2*sum_prev_size+prev_size+t_index);
            *(new_d_ij_buf+2*sum_prev_size+index) = i;
            *(new_d_ij_buf+2*sum_prev_size+new_size+index) = j;
        }
    }
}

__global__ void xor_combinations(
        uint32_t size,
        uint32_t sum_prev_size,
        unsigned char *hashes,
        unsigned char *out,
        uint32_t *ij_buf)
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        uint32_t elem_idx_1 = *(ij_buf+2*sum_prev_size+t_index);
        uint32_t elem_idx_2 = *(ij_buf+2*sum_prev_size+size+t_index);

        // TODO: Skip leading zero bytes
        unsigned char *hash_i = hashes+elem_idx_1*30*sizeof(unsigned char);
        unsigned char *hash_j = hashes+elem_idx_2*30*sizeof(unsigned char);
        unsigned char *hash_out = out+t_index*30*sizeof(unsigned char);

        for(int i=0; i<30; i++) {
            hash_out[i] = hash_i[i] ^ hash_j[i];
        }
    }
}

uint32_t cuda_calculate_pairs(
        unsigned char **d_hashes, 
        uint32_t *d_sort_indices, 
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes) 
{
    uint32_t lookup_size = sizeof(uint32_t)*(1<<24);
    uint32_t hash_count = buf_sizes[r];

    // start lookup table (2**24 elements)
    uint32_t *d_start_lookup;
    cudaMalloc(&d_start_lookup, lookup_size);
    cudaMemset(d_start_lookup, 0, lookup_size); 

    // end lookup table (2**24 elements)
    uint32_t *d_end_lookup;
    cudaMalloc(&d_end_lookup, lookup_size);
    cudaMemset(d_end_lookup, 0, lookup_size); 

    // prefix count lookup table
    uint32_t *d_count_lookup;
    cudaMalloc(&d_count_lookup, lookup_size);  
    cudaMemset(d_count_lookup, 0, lookup_size); 

    // populate start, end, and count arrays
    uint32_t block_count = hash_count / BLOCK_SIZE;
    if(hash_count % BLOCK_SIZE != 0)
        block_count += 1;
    build_start_end_count<<<block_count, BLOCK_SIZE>>>(
        *d_hashes, 
        d_sort_indices, 
        d_start_lookup, 
        d_end_lookup, 
        d_count_lookup, 
        r, 
        hash_count);

    // combination_count lookup table
    uint32_t *d_combination_count_lookup;
    cudaMalloc(&d_combination_count_lookup, lookup_size);  
    cudaMemset(d_combination_count_lookup, 0, lookup_size); 

    // k2
    build_combination_count<<<block_count, BLOCK_SIZE>>>(
        *d_hashes, 
        d_sort_indices, 
        d_start_lookup, 
        d_count_lookup, 
        d_combination_count_lookup, 
        r, 
        hash_count);

    // sum sizes in combination lookup table
    thrust::device_vector<uint32_t> dv_ccl(d_combination_count_lookup, d_combination_count_lookup+(1<<24));

    thrust::inclusive_scan(
        dv_ccl.begin(), dv_ccl.end(), dv_ccl.begin());

    uint32_t *d_combination_sum = thrust::raw_pointer_cast(dv_ccl.data());
    // std::vector<uint32_t> h_combination_sum(1<<24);
    // thrust::copy(dv_ccl.begin(), dv_ccl.end(), h_combination_sum.begin());
    uint32_t *h_combination_sum = (uint32_t *)malloc(1);
    cudaMemcpy(h_combination_sum, d_combination_sum+(1<<24)-1, 1*sizeof(uint32_t), cudaMemcpyDeviceToHost);

    uint32_t size = h_combination_sum[0];
    buf_sizes[r+1] = size;
    // printf("SIZE: %x\n", size);

    printf("ROUND: %x\n", r);
    printf("SIZE: %x\n", size);

    // map prefixes to thread indices
    uint32_t *d_tidx_prefix_map;
    cudaMalloc(&d_tidx_prefix_map, size*sizeof(uint32_t));

    map_index_to_prefix<<<block_count, BLOCK_SIZE>>>(
        *d_hashes, 
        d_sort_indices, 
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

    printf("SIZE_SUM: %x\n", size_sum);

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
    calculate_pairs<<<block_count, BLOCK_SIZE>>>(
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

uint32_t cuda_reduce_pairs(
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes,
        uint32_t size,
        uint32_t sum_prev_size,
        uint32_t block_count) 
{
    // allocate reduction buffer uint32_t * size
    uint32_t *d_reduction;
    cudaMalloc(&d_reduction, size*sizeof(uint32_t));
    // memset to 1
    initialize_to_one<<<block_count, BLOCK_SIZE>>>(d_reduction, size);

    
    uint32_t *d_size_buf;
    cudaMalloc(&d_size_buf, sizeof(uint32_t)*CUDA_K);
    cudaMemset(d_size_buf, 0, sizeof(uint32_t)*CUDA_K);
    cudaMemcpy(d_size_buf, &buf_sizes[1], sizeof(uint32_t)*(r+1), cudaMemcpyHostToDevice);
    
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      fprintf(stderr, "- ERROR: %s \n", cudaGetErrorString(error));
    }

    // uint32_t *reduction = (uint32_t *)malloc(sizeof(uint32_t)*size);
    // cudaMemcpy(reduction, d_reduction, sizeof(uint32_t)*size, cudaMemcpyDeviceToHost);
    // for(int i=0; i<10; i++) {
    //     printf("RED: %x\n", reduction[i]);
    // }

    // clock_t *d_clocks;
    // cudaMalloc(&d_clocks, sizeof(clock_t)*4);

    // write 0 if match
    get_reduction<<<block_count, BLOCK_SIZE>>>(
        //d_clocks,
        d_reduction, 
        r,
        d_size_buf,
        size,
        sum_prev_size,
        *d_ij_buf);


    //clock_t *clocks = (clock_t*)malloc(sizeof(clock_t)*4);
    //cudaMemcpy(clocks, d_clocks, sizeof(clock_t)*4, cudaMemcpyDeviceToHost);
    //for(int i=0; i<4; i++){
    //    printf("C%d: %x\n", i, clocks[i]);
    //}

    // reduce
    thrust::device_vector<uint32_t> dv_reduction(d_reduction, d_reduction+size);
    // suffix_sum reduction buffer
    int reduced_size = thrust::reduce(dv_reduction.begin(), dv_reduction.end());

    if(reduced_size == size)
        return size;

    thrust::exclusive_scan(dv_reduction.begin(), dv_reduction.end(), dv_reduction.begin());
    uint32_t *d_mapping = thrust::raw_pointer_cast(dv_reduction.data());
    
    // uint32_t *h_mapping = (uint32_t *)malloc(dv_reduction.size()*sizeof(uint32_t));
    // cudaMemcpy(h_mapping, d_mapping, dv_reduction.size()*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    // for(int i=0; i<10; i++)
    //     printf("%x\n", h_mapping[i]);

    // uint32_t *h_reduction = (uint32_t *)malloc(sizeof(uint32_t)*size);
    // cudaMemcpy(h_reduction, d_reduction, sizeof(uint32_t)*size, cudaMemcpyDeviceToHost);
    // for(int i=0; i<10; i++)
    //     printf("RED: %x\n", h_reduction[i]);

    // allocate pair buffer, size of sum*2*uint32_t         (sum+1)?
    buf_sizes[r+1] = reduced_size;

    uint32_t size_sum = 0;
    sum_prev_size = 0;
    for(int i=0; i<r+1; i++) {
        size_sum += buf_sizes[1+i];
        if(i != r)
            sum_prev_size += buf_sizes[1+i];
    }

    uint32_t *temp_d_ij_buf;
    cudaMalloc(&temp_d_ij_buf, sizeof(uint32_t)*size_sum*2);
    cudaMemcpy(temp_d_ij_buf, *d_ij_buf, sizeof(uint32_t)*sum_prev_size*2, cudaMemcpyDeviceToDevice);

    remap_reduction<<<block_count, BLOCK_SIZE>>>(
        d_reduction,
        d_mapping,
        *d_ij_buf,
        sum_prev_size,
        size,
        temp_d_ij_buf,
        reduced_size);

    printf("REDUCED_SIZE: %x\n", reduced_size);

    size = reduced_size;
    cudaFree(*d_ij_buf);
    cudaMalloc(d_ij_buf, sizeof(uint32_t)*size_sum*2);
    cudaMemcpy(*d_ij_buf, temp_d_ij_buf, sizeof(uint32_t)*size_sum*2, cudaMemcpyDeviceToDevice);
    cudaFree(temp_d_ij_buf);

    return size;
}

void cuda_xor_step(
        unsigned char **d_hashes, 
        uint32_t *d_sort_indices, 
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes) 
{
    
    uint32_t size = cuda_calculate_pairs(d_hashes, d_sort_indices, d_ij_buf, r, buf_sizes);
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

    if(r != 0) {
        size = cuda_reduce_pairs(d_ij_buf, r, buf_sizes, size, sum_prev_size, block_count);
        size_sum = 0;
        sum_prev_size = 0;
        for(int i=0; i<r+1; i++) {
            size_sum += buf_sizes[1+i];
            if(i != r)
                sum_prev_size += buf_sizes[1+i];
        }

        block_count = size / BLOCK_SIZE;
        if(size % BLOCK_SIZE != 0)
            block_count += 1;
    }

    // xor using pairs and hashes

    // allocate output buf  
    unsigned char *d_output;
    cudaMalloc(&d_output, size*30*sizeof(unsigned char));
    cudaMemset(d_output, 0, size*30*sizeof(unsigned char)); 

    xor_combinations<<<block_count, BLOCK_SIZE>>>(
        size,
        sum_prev_size,
        *d_hashes,
        d_output,
        *d_ij_buf);

    cudaFree(*d_hashes);
    cudaMalloc(d_hashes, size*30*sizeof(unsigned char));
    cudaMemcpy(*d_hashes, d_output, size*30*sizeof(unsigned char), cudaMemcpyDeviceToDevice);
    cudaFree(d_output);


    uint32_t *h_ij_buf = (uint32_t *)malloc(sizeof(uint32_t)*size_sum*2);
    cudaMemcpy(h_ij_buf, *d_ij_buf, sizeof(uint32_t)*size_sum*2, cudaMemcpyDeviceToHost);

    printf("D_IJ_SIZE: %x\n", size_sum*2);

    unsigned char *h_out = (unsigned char *)malloc(size*30*sizeof(unsigned char));
    cudaMemcpy(h_out, *d_hashes, size*30*sizeof(unsigned char), cudaMemcpyDeviceToHost);

    for(int i=0; i<25; i++) {
        unsigned char *hash = h_out+i*30;
        for(int j=0; j<30; j++)
            printf("%02x ", hash[j]);
        uint32_t h_i = *(h_ij_buf+2*sum_prev_size+i);
        uint32_t h_j = *(h_ij_buf+2*sum_prev_size+size+i);
        printf("\t%x, %x", h_i, h_j);
        printf("\n");
    }
}