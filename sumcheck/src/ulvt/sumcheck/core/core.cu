#include <cstdint>
#include <iostream>

#include "../../finite_fields/circuit_generator/unrolled/binary_tower_unrolled.cuh"
#include "../../utils/bitslicing.cuh"
#include "../utils/constants.hpp"
#include "core.cuh"

__host__ __device__ void evaluate_composition_on_batch_row(
	const uint32_t* first_batch_of_row,
	uint32_t* batch_composition_destination,
	const uint32_t composition_size,
	const uint32_t original_evals_per_col
) {
	memcpy(batch_composition_destination, first_batch_of_row, BITS_WIDTH * sizeof(uint32_t));

	for (int operand_in_composition = 1; operand_in_composition < composition_size; ++operand_in_composition) {
		const uint32_t* nth_batch_of_row =
			first_batch_of_row + operand_in_composition * original_evals_per_col * INTS_PER_VALUE;

		multiply_unrolled<TOWER_HEIGHT>(batch_composition_destination, nth_batch_of_row, batch_composition_destination);
	}
}

__host__ __device__ void fold_batch(
	const uint32_t lower_batch[BITS_WIDTH],
	const uint32_t upper_batch[BITS_WIDTH],
	uint32_t dst_batch[BITS_WIDTH],
	const uint32_t coefficient[BITS_WIDTH],
	const bool is_interpolation
) {
	uint32_t xor_of_halves[BITS_WIDTH];

	for (int i = 0; i < BITS_WIDTH; ++i) {
		xor_of_halves[i] = lower_batch[i] ^ upper_batch[i];
	}

	uint32_t product[BITS_WIDTH];
	memset(product, 0, BITS_WIDTH * sizeof(uint32_t));

	// Multiply chunk-wise based on field height of coefficient
	// For random challenges this will be the full 7
	// For interpolation points this will be no more than 2

	if (is_interpolation) {
		for (int i = 0; i < BITS_WIDTH; i += INTERPOLATION_BITS_WIDTH) {
			multiply_unrolled<INTERPOLATION_TOWER_HEIGHT>(xor_of_halves + i, coefficient, product + i);
		}
	} else {
		multiply_unrolled<TOWER_HEIGHT>(xor_of_halves, coefficient, product);
	}

	for (int i = 0; i < BITS_WIDTH; ++i) {
		dst_batch[i] = lower_batch[i] ^ product[i];
	}
}

void fold_small(
	const uint32_t source[BITS_WIDTH],
	uint32_t destination[BITS_WIDTH],
	const uint32_t coefficient[BITS_WIDTH],
	const uint32_t list_len
) {
	uint32_t half_len = list_len / 2;

	uint32_t batch_to_be_multiplied[BITS_WIDTH];

	memcpy(batch_to_be_multiplied, source, BITS_WIDTH * sizeof(uint32_t));

	for (int i = 0; i < BITS_WIDTH; ++i) {
		batch_to_be_multiplied[i] >>= half_len;  // Move the upper half into the lower half of this operand
		batch_to_be_multiplied[i] ^= source[i];  // Add two halves before multiplying
	}

	uint32_t product[BITS_WIDTH];

	multiply_unrolled<TOWER_HEIGHT>(batch_to_be_multiplied, coefficient, product);

	for (int i = 0; i < BITS_WIDTH; ++i) {
		destination[i] = source[i] ^ product[i];
	}
}

template <> 
__global__ void fold_small_kernel<TOWER_HEIGHT>(uint32_t*  d_src,
                                  uint32_t*        d_dst,
                                  uint32_t*  d_coeff,
                                  uint32_t        list_len,
                                  uint32_t        num_cols)
{
    const uint32_t col = blockIdx.x * blockDim.x + threadIdx.x;
    if (col >= num_cols) return;

    const uint32_t half_len = list_len >> 1;

    /* Pointers to the current column -------------------------------- */
    const uint32_t* src = d_src + col * BITS_WIDTH;
    uint32_t*       dst = d_dst + col * BITS_WIDTH;

    /* Local scratch (lives in registers / local memory) ------------- */
    uint32_t batch[BITS_WIDTH];
#pragma unroll
    for (int i = 0; i < BITS_WIDTH; ++i) {
        uint32_t v = src[i];
        batch[i] = (v >> half_len) ^ v;          // bring upper half down + XOR
    }

    uint32_t prod[BITS_WIDTH];
    multiply_unrolled<TOWER_HEIGHT>(batch, d_coeff, prod);

#pragma unroll
    for (int i = 0; i < BITS_WIDTH; ++i)
        dst[i] = src[i] ^ prod[i];               // final xor into destination
}

// __host__ __device__ void compute_sum(
// 	uint32_t sum[INTS_PER_VALUE],
// 	uint32_t bitsliced_batch[BITS_WIDTH],
// 	const uint32_t num_eval_points_being_summed_unpadded
// ) {
// 	BitsliceUtils<BITS_WIDTH>::bitslice_untranspose(bitsliced_batch);

// 	memset(sum, 0, INTS_PER_VALUE * sizeof(uint32_t));

// 	for (uint32_t i = 0; i < min(BITS_WIDTH, INTS_PER_VALUE * num_eval_points_being_summed_unpadded); ++i) {
// 		sum[i % INTS_PER_VALUE] ^= bitsliced_batch[i];
// 	}
// }

// static inline uint32_t parity32(uint32_t x) {
//     x ^= x >> 16;
//     x ^= x >> 8;
//     x ^= x >> 4;
//     x ^= x >> 2;
//     x ^= x >> 1;
//     return x & 1;
// }

// // without unrolling, potentially faster
// __host__ __device__ void compute_sum(
// 	uint32_t sum[INTS_PER_VALUE],
// 	uint32_t bitsliced_batch[BITS_WIDTH],
// 	const uint32_t       num_eval_points
// ) {
//     /* 0. clear the output ------------------------------------------------ */
//     for (uint32_t lane = 0; lane < INTS_PER_VALUE; ++lane)
//         sum[lane] = 0;

//     /* 1. block anatomy --------------------------------------------------- */
//     const uint32_t full_blocks  = num_eval_points >> 5;        /* N / 32 */
//     const uint32_t tail_points  = num_eval_points & 31u;       /* N % 32 */
//     const uint32_t total_blocks = full_blocks + (tail_points ? 1u : 0u);

//     /* 2. iterate over every bit position and every 32-bit limb ---------- */
//     for (uint32_t bit = 0; bit < 32u; ++bit) {
//         for (uint32_t lane = 0; lane < INTS_PER_VALUE; ++lane) {

//             uint32_t bit_parity = 0;

//             /* ---- full 32-point blocks --------------------------------- */
//             for (uint32_t blk = 0; blk < full_blocks; ++blk) {

// #if BIT_MAJOR_LAYOUT
//                 uint32_t slice =
//                     bitsliced_batch[ blk*BITS_WIDTH + bit*INTS_PER_VALUE + lane ];
// #else
//                 uint32_t slice =
//                     bitsliced_batch[ blk*BITS_WIDTH + lane*32u + bit ];
// #endif
//                 bit_parity ^= parity32(slice);                /* XOR parity */
//             }

//             /* ---- tail block (0 < tail_points < 32) -------------------- */
//             if (tail_points) {

// #if BIT_MAJOR_LAYOUT
//                 uint32_t slice =
//                     bitsliced_batch[ full_blocks*BITS_WIDTH + bit*INTS_PER_VALUE + lane ];
// #else
//                 uint32_t slice =
//                     bitsliced_batch[ full_blocks*BITS_WIDTH + lane*32u + bit ];
// #endif
//                 const uint32_t mask = (1u << tail_points) - 1u;
//                 bit_parity ^= parity32(slice & mask);
//             }

//             /* ---- drop the single parity bit into the result limb ------ */
//             sum[lane] ^= bit_parity << bit;                   /* XOR – safe for reuse */
//         }
//     }
// }

/*****************************************************************************
*  compute_sum.cu  –  CUDA version that uses the native POPC unit            *
*                                                                            *
*  • Handles any value width:   value_bits = 32 × INTS_PER_VALUE             *
*  • Handles any batch size N   (works for N ≤ 32 and N > 32)                *
*  • No atomics, no shared-memory reduction                                  *
*  • Parity is computed with a single POPC instruction per slice             *
*                                                                            *
*  Pick **exactly one** layout flag below to match how your slices are laid   *
*  out in memory.                                                             *
*****************************************************************************/

#define INTS_PER_VALUE    4      // 128-bit value  → 4 × 32-bit limbs
#define BIT_MAJOR_LAYOUT  0      // slice index =  block*W + bit*L + lane
#define LANE_MAJOR_LAYOUT 1      // slice index =  block*W + lane*32 + bit
static_assert(BIT_MAJOR_LAYOUT ^ LANE_MAJOR_LAYOUT, "pick exactly one");

#define BITS_PER_LIMB  32
#define BITS_WIDTH    (BITS_PER_LIMB * INTS_PER_VALUE)

/* --------------------------------------------------------------------- */
/* Kernel: one thread per bit, one block per limb                        */
/* --------------------------------------------------------------------- */
__global__ void compute_sum_kernel(uint32_t      *sum_out,   // OUT: INTS_PER_VALUE words
                                   const uint32_t *slices,   // IN : ≥ ceil(N/32) * BITS_WIDTH words
                                   uint32_t        N)        // IN : evaluation-point count (≥1)
{
    const uint32_t lane   = blockIdx.x;    // which 32-bit limb (0 … INTS_PER_VALUE-1)
    const uint32_t bit    = threadIdx.x;   // which bit position (0 … 31)
    const uint32_t blocks = (N + 31) >> 5; // ceil(N/32)
    const uint32_t tail   =  N & 31u;      // N mod 32   (0 → no tail block)

    uint32_t parity = 0;

    /* ---- all full 32-point blocks ----------------------------------- */
    for (uint32_t blk = 0; blk < blocks - (tail ? 1u : 0u); ++blk) {
#if BIT_MAJOR_LAYOUT
        uint32_t slice = slices[ blk*BITS_WIDTH + bit*INTS_PER_VALUE + lane ];
#else
        uint32_t slice = slices[ blk*BITS_WIDTH + lane*32u + bit ];
#endif
        parity ^= (__popc(slice) & 1u);             // POPC gives 32-bit parity in 1 inst
    }

    /* ---- optional tail block (mask out padding bits) ---------------- */
    if (tail) {
#if BIT_MAJOR_LAYOUT
        uint32_t slice = slices[ (blocks-1u)*BITS_WIDTH + bit*INTS_PER_VALUE + lane ];
#else
        uint32_t slice = slices[ (blocks-1u)*BITS_WIDTH + lane*32u + bit ];
#endif
        slice &= (1u << tail) - 1u;                 // keep only the real N % 32 bits
        parity ^= (__popc(slice) & 1u);
    }

    /* ---- warp-wide pack: thread-bit → result word ------------------- */
    uint32_t limb_word = __ballot_sync(0xFFFFFFFFu, parity);

    /* ---- one thread per block writes its limb ----------------------- */
    if (threadIdx.x == 0)
        sum_out[lane] = limb_word;                  // no atomics – block owns this limb
}

/* --------------------------------------------------------------------- */
/* Convenience launcher (host side)                                      */
/* --------------------------------------------------------------------- */
void compute_sum(
	uint32_t sum[INTS_PER_VALUE],
	uint32_t bitsliced_batch[BITS_WIDTH],
	const uint32_t       num_eval_points)     /* N ≥ 1 */
{
    /* 1. Work geometry ------------------------------------------------- */
    const uint32_t blocks_host = (num_eval_points + 31) >> 5;   /* ceil(N/32) */
    const size_t   slice_words = static_cast<size_t>(blocks_host) * BITS_WIDTH;

    /* 2. Allocate device memory --------------------------------------- */
    uint32_t *d_slices = nullptr, *d_sum = nullptr;
    cudaMalloc(&d_slices, slice_words * sizeof(uint32_t));
    cudaMalloc(&d_sum,    INTS_PER_VALUE * sizeof(uint32_t));

    /* 3. Copy slices to the device ------------------------------------ */
    cudaMemcpy(d_slices, bitsliced_batch,
               slice_words * sizeof(uint32_t), cudaMemcpyHostToDevice);

    /* 4. Launch: one 32-thread block per limb ------------------------- */
    dim3 grid (INTS_PER_VALUE);
    dim3 block(32);
    compute_sum_kernel<<<grid, block>>>(d_sum, d_slices, num_eval_points);

    /* 5. Copy the result back ----------------------------------------- */
    cudaMemcpy(sum, d_sum,
               INTS_PER_VALUE * sizeof(uint32_t), cudaMemcpyDeviceToHost);

    /* 6. Cleanup ------------------------------------------------------- */
    cudaFree(d_slices);
    cudaFree(d_sum);
}