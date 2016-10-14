// Copyright Douglas Goddard 2016
// Licensed under the MIT license

#include "csolv.h"
#include <stdio.h>

#define IV_0 0x6a09e667f3bcc908UL
#define IV_1 0xbb67ae8584caa73bULL
#define IV_2 0x3c6ef372fe94f82bULL
#define IV_3 0xa54ff53a5f1d36f1ULL
#define IV_4 0x510e527fade682d1ULL
#define IV_5 0x9b05688c2b3e6c1fULL
#define IV_6 0x1f83d9abfb41bd6bULL
#define IV_7 0x5be0cd19137e2179ULL

__device__ static const uint8_t blake2b_sigma[12][16] = {
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 } ,
    { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 } ,
    { 11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4 } ,
    {  7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8 } ,
    {  9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13 } ,
    {  2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9 } ,
    { 12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11 } ,
    { 13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10 } ,
    {  6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5 } ,
    { 10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13 , 0 } ,
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 } ,
    { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 }
};

// __device__ static inline uint64_t rotr64(const uint64_t w, const unsigned c) {
//     return (w >> c) | (w << (64-c));
// }

__device__ static inline void G1(uint8_t r, uint8_t i, uint64_t *m, uint64_t *v) {
    asm("{ \n\t"
        "  .reg .u64 t1; \n\t"
        "  add.u64 %0, %0, %1; \n\t"
        "  add.u64 %0, %0, %4; \n\t"
        "  xor.b64 %3, %3, %0; \n\t"
        "  shr.b64 t1, %3, 32; \n\t"
        "  shl.b64 %3, %3, 32; \n\t"
        "  or.b64  %3, t1, %3; \n\t"
        "  add.u64 %2, %2, %3; \n\t"
        "  xor.b64 %1, %1, %2; \n\t"
        "  shr.b64 t1, %1, 24; \n\t"
        "  shl.b64 %1, %1, 40; \n\t"
        "  or.b64  %1, t1, %1; \n\t"
        "  add.u64 %0, %0, %1; \n\t"
        "  add.u64 %0, %0, %5; \n\t"
        "  xor.b64 %3, %3, %0; \n\t"
        "  shr.b64 t1, %3, 16; \n\t"
        "  shl.b64 %3, %3, 48; \n\t"
        "  or.b64  %3, t1, %3; \n\t"
        "  add.u64 %2, %2, %3; \n\t"
        "  xor.b64 %1, %1, %2; \n\t"
        "  shr.b64 t1, %1, 63; \n\t"
        "  shl.b64 %1, %1,  1; \n\t"
        "  or.b64  %1, t1, %1; \n\t"
        "} "
        : "+l"(v[0 +i]),
        "+l"(v[4 +i]),
        "+l"(v[8 +i]),
        "+l"(v[12+i]),
        "+l"(m[blake2b_sigma[r][2*i+0]]),
        "+l"(m[blake2b_sigma[r][2*i+1]]) :);
}

__device__ static inline void G2(uint8_t r, uint8_t i, uint64_t *m, uint64_t *v) {
    asm("{ \n\t"
        "  .reg .u64 t1; \n\t"
        "  add.u64 %0, %0, %1; \n\t"
        "  add.u64 %0, %0, %4; \n\t"
        "  xor.b64 %3, %3, %0; \n\t"
        "  shr.b64 t1, %3, 32; \n\t"
        "  shl.b64 %3, %3, 32; \n\t"
        "  or.b64  %3, t1, %3; \n\t"
        "  add.u64 %2, %2, %3; \n\t"
        "  xor.b64 %1, %1, %2; \n\t"
        "  shr.b64 t1, %1, 24; \n\t"
        "  shl.b64 %1, %1, 40; \n\t"
        "  or.b64  %1, t1, %1; \n\t"
        "  add.u64 %0, %0, %1; \n\t"
        "  add.u64 %0, %0, %5; \n\t"
        "  xor.b64 %3, %3, %0; \n\t"
        "  shr.b64 t1, %3, 16; \n\t"
        "  shl.b64 %3, %3, 48; \n\t"
        "  or.b64  %3, t1, %3; \n\t"
        "  add.u64 %2, %2, %3; \n\t"
        "  xor.b64 %1, %1, %2; \n\t"
        "  shr.b64 t1, %1, 63; \n\t"
        "  shl.b64 %1, %1,  1; \n\t"
        "  or.b64  %1, t1, %1; \n\t"
        "} "
        : "+l"(v[0 +(i-4)]),
        "+l"(v[4 +((i-3)%4)]),
        "+l"(v[8 +((i-2)%4)]),
        "+l"(v[12+((i-1)%4)]),
        "+l"(m[blake2b_sigma[r][2*i+0]]),
        "+l"(m[blake2b_sigma[r][2*i+1]]) :);
}

__device__ void blake2b_compress(
        uint64_t *s_h, 
        uint64_t s_t, 
        uint64_t s_f, 
        const uint8_t block[CUDA_BLAKE2B_BLOCKBYTES] ) 
{
    uint64_t m[16];
    uint64_t v[16];

    #define SET_M(i) \
        do { \
            m[i] = *(uint64_t *)(block + i * sizeof(uint64_t)); \
        } while(0)

    SET_M(0);
    SET_M(1);
    SET_M(2);
    SET_M(3);
    SET_M(4);
    SET_M(5);
    SET_M(6);
    SET_M(7);
    SET_M(8);
    SET_M(9);
    SET_M(10);
    SET_M(11);
    SET_M(12);
    SET_M(13);
    SET_M(14);
    SET_M(15);

    #define SET_V(i) \
        do { \
            v[i] = s_h[i]; \
        } while(0)

    SET_V(0);
    SET_V(1);
    SET_V(2);
    SET_V(3);
    SET_V(4);
    SET_V(5);
    SET_V(6);
    SET_V(7);

    v[ 8] = IV_0;
    v[ 9] = IV_1;
    v[10] = IV_2;
    v[11] = IV_3;
    v[12] = s_t ^ IV_4;
    v[13] = 0 ^ IV_5;
    v[14] = s_f ^ IV_6;
    v[15] = 0 ^ IV_7;

    #define ROUND(r)  \
    do { \
        G1(r,0,&m[0],&v[0]); \
        G1(r,1,&m[0],&v[0]); \
        G1(r,2,&m[0],&v[0]); \
        G1(r,3,&m[0],&v[0]); \
        G2(r,4,&m[0],&v[0]); \
        G2(r,5,&m[0],&v[0]); \
        G2(r,6,&m[0],&v[0]); \
        G2(r,7,&m[0],&v[0]); \
    } while(0)

    ROUND(0);
    ROUND(1);
    ROUND(2);
    ROUND(3);
    ROUND(4);
    ROUND(5);
    ROUND(6);
    ROUND(7);
    ROUND(8);
    ROUND(9);
    ROUND(10);
    ROUND(11);

    #define SET_H(i) \
        do { \
            s_h[i] = s_h[i] ^ v[i] ^ v[i + 8]; \
        } while(0);

    SET_H(0);
    SET_H(1);
    SET_H(2);
    SET_H(3);
    SET_H(4);
    SET_H(5);
    SET_H(6);
    SET_H(7);

    #undef SET_M
    #undef SET_V
    #undef SET_H
    #undef ROUND
}

__global__ void hash_many(
        uint8_t *d_s_buf, 
        size_t s_buflen,
        unsigned char *d_out) 
{
    uint32_t offset = blockDim.x * blockIdx.x + threadIdx.x;
    // this is just the IV, everything else is in buf
    uint64_t s_h[8] = { 0x6a09e667f2bdc93a ,  0xbb67ae8584caa73bL ,  0x3c6ef372fe94f82b ,  0xa54ff53a5f1d36f1L ,  0x510e527fade682d1 ,  0x9b05688c2b3e6c1fL ,  0x48ec89c3fb41bd62 ,  0x5be0cd19137e2179 };
    uint64_t s_t = 0x0;
    uint64_t s_f = 0x0;
    uint8_t s_buf[256];

    memcpy(&s_buf[0], d_s_buf, 128*2*sizeof(uint8_t));

    // update 
    uint32_t le_offset = htole32(offset);
    memcpy(s_buf+s_buflen, &le_offset, sizeof(uint32_t));
    s_buflen += 4;

    // final
    s_t += CUDA_BLAKE2B_BLOCKBYTES;
    blake2b_compress(s_h, s_t, s_f, s_buf);
    s_buflen -= CUDA_BLAKE2B_BLOCKBYTES;
    memcpy(s_buf, s_buf+CUDA_BLAKE2B_BLOCKBYTES, s_buflen);

    s_t += s_buflen;
    s_f = (uint64_t)-1;
    memset(s_buf+s_buflen, 0, 2*CUDA_BLAKE2B_BLOCKBYTES-(s_buflen));
    blake2b_compress(s_h, s_t, s_f, s_buf);

    uint32_t hash_len = (512/CUDA_N)*CUDA_N/8;
    memcpy(d_out+offset*hash_len, &s_h[0], hash_len);
}

void cuda_generate_hash(
        const CUDA_BLAKE2B_STATE *base_state, 
        uint32_t hash_count, 
        unsigned char *h_out)
{
    if(hash_count % BLOCK_SIZE != 0) {
        printf("error: Invalid hash count!\n");
        exit(1);
    }

    uint32_t hash_len = (512/CUDA_N)*CUDA_N/8;
    uint32_t out_size = hash_len * hash_count * sizeof(unsigned char);

    uint8_t *d_s_buf;
    unsigned char *d_out;
    
    cudaMalloc(&d_s_buf, sizeof(uint8_t)*2*128);
    cudaMalloc(&d_out, out_size);
    
    cudaMemcpy(d_s_buf, base_state->buf, sizeof(uint8_t)*128*2, cudaMemcpyHostToDevice);

    uint32_t block_count = hash_count / BLOCK_SIZE;
    hash_many<<<block_count, BLOCK_SIZE>>>(d_s_buf, base_state->buflen, d_out);

    // wait until tasks are completed
    cudaDeviceSynchronize();

    // check for errors
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    }

    cudaMemcpy(h_out, d_out, out_size, cudaMemcpyDeviceToHost);

    cudaFree(d_s_buf);
    cudaFree(d_out);
}
