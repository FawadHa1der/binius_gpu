#ifndef LINEAR_CHECK_H
#define LINEAR_CHECK_H
#include "field.h"
#include "types.h"
#include "mle_poly.h"
#include "compressed_poly.h"
#include "evaluations.h"
#include "prod_check.h"
#include "univariate_poly.h"
#include "matrix_linear.h"





typedef struct {
    const MatrixLinear *matrix;       // Linear transformation matrix
    const MLE_POLY_SEQUENCE *polys;           // Array of pointers to input polynomials
    const Points *points;             // Points for polynomial evaluation
    size_t num_vars;                  // Total number of variables
    size_t num_active_vars;           // Number of active variables
    F128 *initial_claims;             // Array of M initial claims
    size_t N;                         // Number of input polynomials
    size_t M;                         // Number of output claims
} LinCheckBuilder;

LinCheckBuilder* lincheck_builder_new(
    MLE_POLY_SEQUENCE *polys,
    const Points *points,
    const MatrixLinear *matrix,
    size_t num_active_vars,
    const F128 *initial_claims,
    size_t N,
    size_t M
);

ProdCheck* lincheck_builder_build(LinCheckBuilder *builder, const F128 *gamma);
void lincheck_builder_free(LinCheckBuilder *builder);

#endif