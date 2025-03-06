#ifndef MATRX_UTILS_H
#define MATRX_UTILS_H


#include <stdint.h>
#include <stdbool.h>
#include <stddef.h> // for size_t
#include  "field.h"

typedef struct {
    F128 cols[128];
} Matrix;

/* Functions that mimic the Rust 'Matrix' impl */
Matrix matrix_new(const F128 *cols /* length=128 */);
F128 matrix_apply(const Matrix *m, F128 vec);
Matrix matrix_compose(const Matrix *m, const Matrix *b);
Matrix matrix_diag(void);
void matrix_swap_cols(Matrix *m, size_t i, size_t j);
void matrix_triang(Matrix *m, size_t i, size_t j);

/* 
 * matrix_inverse => returns a bool for success (true) or failure (false). 
 * If success, *out_inv is the inverse.
 */
bool matrix_inverse(const Matrix *in_m, Matrix *out_inv);

/*
 * 3) Additional free functions
 */
unsigned log2_exact(unsigned long long x);

// Like "u128_idx(x, i)", but for a F128's bits
bool f128_idx(const F128 *x, unsigned i);

// Convert a F128 into an array of 128 bool bits (little-endian).
void binary128_to_bits(F128 x, bool out_bits[128]);

// Reconstruct a 128-bit integer from an array of bits
void u128_from_bits(const bool *bits, size_t length, uint64_t *low, uint64_t *high);

void f128_get(const F128 *x, uint64_t *low, uint64_t *high);

#endif // MATRX_UTILS_H