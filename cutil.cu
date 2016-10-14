// Copyright Douglas Goddard 2016
// Licensed under the MIT license

#include "csolv.h"

// this is specific to N=200, K=9
__global__ void expand_array(
        unsigned char *d_in, 
        unsigned char *d_out) 
{
    uint32_t offset = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned char *input = d_in+offset*5*sizeof(unsigned char);
    unsigned char *output = d_out+offset*6*sizeof(unsigned char);

    output[0] = input[0] >> 4;
    output[1] = input[0] << 4 | input[1] >> 4;
    output[2] = input[1] << 4 | input[2] >> 4;
    output[3] = input[2] & 0xf;
    output[4] = input[3];
    output[5] = input[4];
}

__global__ void get_sort_keys(
        unsigned char* d_hashes, 
        uint32_t *d_sort_keys, 
        uint32_t *d_sort_indices,
        uint32_t r,
        uint32_t num_keys) 
{
    uint32_t index = blockDim.x * blockIdx.x + threadIdx.x;
    if(index < num_keys) {
        d_sort_indices[index] = index;

        unsigned char* input = d_hashes+index*30*sizeof(unsigned char)+3*r;
        uint32_t sort_key = input[0] << 16 | input[1] << 8 | input[2];
        d_sort_keys[index] = sort_key;
    }
}

void cuda_sort_step(
        unsigned char* d_hashes, 
        uint32_t *d_sort_indices_out, 
        uint32_t r, 
        uint32_t num_keys) 
{

    uint32_t *d_sort_keys;
    uint32_t *d_sort_indices;

    cudaMalloc(&d_sort_keys, sizeof(uint32_t)*num_keys);
    cudaMalloc(&d_sort_indices, sizeof(uint32_t)*num_keys);

    uint32_t block_count = num_keys / BLOCK_SIZE;
    if(num_keys % BLOCK_SIZE != 0)
        block_count += 1;
    get_sort_keys<<<block_count, BLOCK_SIZE>>>(d_hashes, d_sort_keys, d_sort_indices, r, num_keys);

    thrust::device_vector<uint32_t> dv_sort_keys(d_sort_keys, d_sort_keys+num_keys);
    thrust::device_vector<uint32_t> dv_sort_indices(d_sort_indices, d_sort_indices+num_keys);
    
    cudaFree(d_sort_keys);
    cudaFree(d_sort_indices);
    
    thrust::sort_by_key(dv_sort_keys.begin(), dv_sort_keys.end(), dv_sort_indices.begin());

    thrust::copy(dv_sort_indices.begin(), dv_sort_indices.end(), d_sort_indices_out);
}

__global__ void initialize_to_one(
        uint32_t *reduction,
        uint32_t size) 
{
    uint32_t t_index = blockDim.x * blockIdx.x + threadIdx.x;
    if(t_index < size) {
        reduction[t_index] = 1;
    }
}

void recover_indices_recursive(
    std::vector<uint32_t> *indices, 
    uint32_t **i_bufs, 
    uint32_t **j_bufs, 
    uint32_t r, 
    uint32_t index) 
{
    if(r == 0) {
        printf("ADDING:\t%x, %x\n", i_bufs[r][index], j_bufs[r][index]);
        indices->push_back(i_bufs[r][index]);
        indices->push_back(j_bufs[r][index]);
    } else {
        uint32_t i = i_bufs[r][index];
        uint32_t j = j_bufs[r][index];
        printf("index:\t%x\n", index);
        printf("i, j:\t%x, %x\n", i, j);
        recover_indices_recursive(indices, i_bufs, j_bufs, r-1, i);
        recover_indices_recursive(indices, i_bufs, j_bufs, r-1, j);
    }
}

void cuda_expand_array(
        unsigned char* h_in, 
        uint32_t in_size,
        unsigned char **d_out,
        uint32_t out_size)
{
    unsigned char *d_in;
    cudaMalloc(&d_in, in_size);
    
    cudaMemcpy(d_in, h_in, sizeof(unsigned char)*in_size, cudaMemcpyHostToDevice);

    uint32_t thread_count = in_size / 5;
    uint32_t block_count = thread_count / BLOCK_SIZE;
    expand_array<<<block_count, BLOCK_SIZE>>>(d_in, *d_out);
    cudaFree(d_in);
}

void cuda_solve_hashes(
        unsigned char* h_in, 
        uint32_t in_size,
        unsigned char *h_out,
        uint32_t out_size) 
{
    unsigned char *d_out;
    cudaMalloc(&d_out, out_size);

    cuda_expand_array(h_in, in_size, &d_out, out_size);

    // copy out hashes so we can reuse d_out in loop
    cudaMemcpy(h_out, d_out, out_size, cudaMemcpyDeviceToHost);

    // track sizes
    uint32_t *buf_sizes = (uint32_t *)malloc(sizeof(uint32_t)*(CUDA_K+1));
    buf_sizes[0] = (1<<21);

    uint32_t *d_ij_buf;
    uint32_t *d_sort_indices; 
    // for(int r=0; r<3; r++) {
    for(int r=0; r<CUDA_K; r++) {
        // realloc sort_indices to size
        cudaMalloc(&d_sort_indices, sizeof(uint32_t)*buf_sizes[r]);
        cuda_sort_step(d_out, d_sort_indices, r, buf_sizes[r]);
        if(r < (CUDA_K-1))
            cuda_xor_step(&d_out, d_sort_indices, &d_ij_buf, r, buf_sizes);
        else 
           cuda_final_collision_step(&d_out, d_sort_indices, &d_ij_buf, r, buf_sizes);
        cudaFree(d_sort_indices);
    }
    
    // std::vector<uint32_t> indices;
    // recover_indices_recursive(&indices, i_bufs, j_bufs, 3, 0x1bf94);

    // for(int i=0; i<indices.size(); i++) {
    //     printf("%x\n", indices[i]);
    // }


    // uint32_t *h_sort_indices = (uint32_t *)malloc(sizeof(uint32_t)*(1<<21));
    // cudaMemcpy(h_sort_indices, d_sort_indices, sizeof(uint32_t)*(1<<21), cudaMemcpyDeviceToHost);

    // printf("HI: %x\n", h_sort_indices[0]);
    // printf("HI: %x\n", h_sort_indices[1]);
    // printf("HI: %x\n", h_sort_indices[2]);

    
    // for(int i=0; i<25; i++) {
    //     unsigned char *hash = h_out+h_sort_indices[i]*30;
    //     for(int j=0; j<30; j++)
    //         printf("%02x ", hash[j]);
    //     printf("\n");
    // }


    cudaFree(d_out);
}
