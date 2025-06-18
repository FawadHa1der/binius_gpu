#pragma once
#include <cstdint>

#include "../utils/constants.hpp"

__host__ __device__ void evaluate_composition_on_batch_row(
	const uint32_t* first_batch_of_row,
	uint32_t* batch_composition_destination,
	const uint32_t composition_size,
	const uint32_t original_evals_per_col
);

__host__ __device__ void fold_batch(
	const uint32_t lower_batch[BITS_WIDTH],
	const uint32_t upper_batch[BITS_WIDTH],
	uint32_t dst_batch[BITS_WIDTH],
	const uint32_t coefficient[BITS_WIDTH],
	const bool is_interpolation
);

void fold_small(
	const uint32_t source[BITS_WIDTH],
	uint32_t destination[BITS_WIDTH],
	const uint32_t coefficient[BITS_WIDTH],
	const uint32_t list_len
);

__host__ __device__ void compute_sum(
	uint32_t sum[INTS_PER_VALUE],
	uint32_t bitsliced_batch[BITS_WIDTH],
	const uint32_t num_eval_points_being_summed_unpadded
);

template <int HEIGHT>
__global__ void fold_small_kernel(uint32_t*  d_src,
                                  uint32_t*        d_dst,
                                  uint32_t*  d_coeff,
                                  uint32_t        list_len,
                                  uint32_t        num_cols);

__host__ __device__ void compute_sum_gpu(
    uint32_t sum[INTS_PER_VALUE],
    uint32_t bitsliced_batch[BITS_WIDTH],
    const uint32_t num_eval_points_being_summed_unpadded
);

template<int HEIGHT>
__global__ void evaluate_composition_kernel(const uint32_t* __restrict__ d_src,
                                            uint32_t*       __restrict__ d_dst,
                                            uint32_t        composition_size,
                                            uint32_t        original_evals_per_col,
                                            uint32_t        row_stride,        // words
                                            uint32_t        num_rows);

void evaluate_composition_gpu(const uint32_t* h_batches,   // all rows, on host
                              uint32_t*       h_results,   // one row → BITS_WIDTH words
                              uint32_t        composition_size,
                              uint32_t        original_evals_per_col,
                              uint32_t        num_rows);

