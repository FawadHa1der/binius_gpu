#ifndef EFFICIENT_MATRIX_H
#define EFFICIENT_MATRIX_H

#include "types.h"
#include "field.h"
#include "mle_poly.h"


#define NUM_ROWS 128
#define NUM_COLS 128

// We'll define that 256*16 => 4096
#define EFFICIENT_MATRIX_SIZE (256 * 16)

// The "EfficientMatrix" is an array of 4096 F128
typedef struct {
    F128 data[EFFICIENT_MATRIX_SIZE];
} EfficientMatrix;


EfficientMatrix* from_rows(const F128 rows[NUM_ROWS]);
EfficientMatrix* from_cols(const F128 cols[NUM_COLS]);
EfficientMatrix* from_frobenius_inv_lc(const F128 gammas[NUM_COLS]);
EfficientMatrix* efficient_matrix_zero(void);
F128 efficient_matrix_apply(const EfficientMatrix *matrix, F128 rhs);

void efficient_matrix_free(EfficientMatrix *matrix);


#endif