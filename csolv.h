#include <sodium.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <set>

#define BLOCK_SIZE 128
#define CUDA_BLAKE2B_BLOCKBYTES 128
#define CUDA_N 200
#define CUDA_K 9
#define CUDA_BLAKE2B_STATE crypto_generichash_blake2b_state

void cuda_generate_hash(
    const CUDA_BLAKE2B_STATE *base_state, 
    uint32_t hash_count, 
    unsigned char *h_out);

void cuda_solve_hashes(
        unsigned char* h_in, 
        uint32_t in_size,
        unsigned char *h_out,
        uint32_t out_size);

void cuda_final_collision_step(
        unsigned char **d_hashes, 
        uint32_t *d_sort_indices, 
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes);

void cuda_xor_step(
        unsigned char **d_hashes, 
        uint32_t *d_sort_indices, 
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes);

uint32_t cuda_reduce_pairs(
        uint32_t **d_ij_buf,
        uint32_t r, 
        uint32_t *buf_sizes,
        uint32_t size,
        uint32_t sum_prev_size,
        uint32_t block_count);

__global__ void initialize_to_one(
        uint32_t *reduction,
        uint32_t size);