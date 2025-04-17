#include "test_linear_check.h"
#include "unity.h"
#include "../field.h"
#include "../matrix_utils.h"
#include "../matrix_linear.h"
#include "../prod_check.h"
#include "../compressed_poly.h"
#include "../univariate_poly.h"
#include "../linear_check.h"
#include "../debug_utils.h"


void test_default_lincheck() {
    Points points = points_default(); // assume len = 0
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 1); // N=2, len=1 each
    MatrixLinear matrix = matrix_linear_new(0, 0, NULL, 0);
    F128 initial_claims[2] = { f128_zero(), f128_zero() };

    LinCheckBuilder* builder = lincheck_builder_new(polys, &points, &matrix, 0, initial_claims, 2, 2);

    assert(builder->num_vars == 0);
    assert(builder->num_active_vars == 0);
    assert(builder->matrix->n_in == 0);
    assert(builder->matrix->n_out == 0);
    for (size_t i = 0; i < 2; i++) {
        assert(f128_eq(builder->initial_claims[i], f128_zero()));
    }

    lincheck_builder_free(builder);
    mle_sequence_free(polys);
}