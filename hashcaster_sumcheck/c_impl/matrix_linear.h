#ifndef MATRIX_LINEAR_H
#define MATRIX_LINEAR_H
#include <stdint.h>
#include <stddef.h> // for size_t
#include "field.h"
#include "matrix_utils.h"

typedef struct {
    size_t n_in;   // number of input dimensions (columns)
    size_t n_out;  // number of output dimensions (rows)
    F128* entries; // flattened matrix of length (n_in * n_out)
} MatrixLinear;


// A "constructor" that returns a new MatrixLinear
// The length of entries must be n_in*n_out.
MatrixLinear matrix_linear_new(size_t n_in, size_t n_out, const F128* entries, size_t len_entries);

// Freed
void matrix_linear_free(MatrixLinear *mat);

// Return the dimension
size_t matrix_linear_n_in(const MatrixLinear *mat);
size_t matrix_linear_n_out(const MatrixLinear *mat);

// apply: mat × input => output
// input has length = mat->n_in
// output has length= mat->n_out
void matrix_linear_apply(const MatrixLinear *mat,
                         const F128 *input, size_t in_len,
                         F128 *output, size_t out_len);

// apply_transposed: mat^T × input => output
// input has length= mat->n_out
// output has length= mat->n_in
void matrix_linear_apply_transposed(const MatrixLinear *mat,
                                    const F128 *input, size_t in_len,
                                    F128 *output, size_t out_len);

#endif