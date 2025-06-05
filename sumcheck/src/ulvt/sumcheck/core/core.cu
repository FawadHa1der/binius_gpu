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

static inline uint32_t parity32(uint32_t x) {
    x ^= x >> 16;
    x ^= x >> 8;
    x ^= x >> 4;
    x ^= x >> 2;
    x ^= x >> 1;
    return x & 1;
}

__host__ __device__ void compute_sum(
	uint32_t sum[INTS_PER_VALUE],
	uint32_t bitsliced_batch[BITS_WIDTH],
	const uint32_t       num_eval_points
) {
    /* 0. clear the output ------------------------------------------------ */
    for (uint32_t lane = 0; lane < INTS_PER_VALUE; ++lane)
        sum[lane] = 0;

    /* 1. block anatomy --------------------------------------------------- */
    const uint32_t full_blocks  = num_eval_points >> 5;        /* N / 32 */
    const uint32_t tail_points  = num_eval_points & 31u;       /* N % 32 */
    const uint32_t total_blocks = full_blocks + (tail_points ? 1u : 0u);

    /* 2. iterate over every bit position and every 32-bit limb ---------- */
    for (uint32_t bit = 0; bit < 32u; ++bit) {
        for (uint32_t lane = 0; lane < INTS_PER_VALUE; ++lane) {

            uint32_t bit_parity = 0;

            /* ---- full 32-point blocks --------------------------------- */
            for (uint32_t blk = 0; blk < full_blocks; ++blk) {

#if BIT_MAJOR_LAYOUT
                uint32_t slice =
                    bitsliced_batch[ blk*BITS_WIDTH + bit*INTS_PER_VALUE + lane ];
#else
                uint32_t slice =
                    bitsliced_batch[ blk*BITS_WIDTH + lane*32u + bit ];
#endif
                bit_parity ^= parity32(slice);                /* XOR parity */
            }

            /* ---- tail block (0 < tail_points < 32) -------------------- */
            if (tail_points) {

#if BIT_MAJOR_LAYOUT
                uint32_t slice =
                    bitsliced_batch[ full_blocks*BITS_WIDTH + bit*INTS_PER_VALUE + lane ];
#else
                uint32_t slice =
                    bitsliced_batch[ full_blocks*BITS_WIDTH + lane*32u + bit ];
#endif
                const uint32_t mask = (1u << tail_points) - 1u;
                bit_parity ^= parity32(slice & mask);
            }

            /* ---- drop the single parity bit into the result limb ------ */
            sum[lane] ^= bit_parity << bit;                   /* XOR – safe for reuse */
        }
    }
}