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
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 1, f128_zero()); // N=2, len=1 each
    MatrixLinear* matrix = matrix_linear_new(0, 0, NULL, 0);
    F128 initial_claims[2] = { f128_zero(), f128_zero() };

    LinCheckBuilder* builder = lincheck_builder_new(polys, &points, matrix, 0, initial_claims, 2, 2);

    assert(builder->num_vars == 0);
    assert(builder->num_active_vars == 0);
    assert(builder->matrix->n_in == 0);
    assert(builder->matrix->n_out == 0);
    for (size_t i = 0; i < 2; i++) {
        assert(f128_eq(builder->initial_claims[i], f128_zero()));
    }

    lincheck_builder_free(builder);
    mle_sequence_free(polys);
    matrix_linear_free(matrix);
}


void test_new_lincheck() {
    size_t num_vars = 3;
    Points* points = points_init(num_vars, f128_zero());
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(1, 1 << num_vars, f128_zero());

    Points* matrix_points = points_init(16, f128_zero());
    MatrixLinear* matrix = matrix_linear_new(4, 4, matrix_points->elems, matrix_points->len);
    F128 initial_claims[1] = { f128_zero() };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, 2, initial_claims, 1, 1);
    assert(builder->num_vars == points->len);
    assert(builder->num_active_vars == 2);
    assert(matrix->n_in == 4);
    assert(matrix->n_out == 4);

    lincheck_builder_free(builder);
    mle_sequence_free(polys);
}


void test_invalid_matrix_dimensions() {
    Points *points = points_init(0, f128_zero());

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(1, 1, f128_zero());
    Points* matrix_points = points_init(16, f128_zero());
    MatrixLinear* matrix = matrix_linear_new(3, 4, matrix_points->elems, matrix_points->len);
    F128 initial_claims[1] = { f128_zero() };

    // This should assert fail
    lincheck_builder_new(polys, points, matrix, 2, initial_claims, 1, 1);
}


void test_invalid_polynomial_length() {
    Points *points = points_init(0, f128_zero());

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(1, 1, f128_zero()); // invalid length for 3 active vars
    Points* matrix_points = points_init(16, f128_zero());
    MatrixLinear *matrix = matrix_linear_new(4, 4, matrix_points->elems, matrix_points->len);
    F128 initial_claims[1] = { f128_zero() };

    lincheck_builder_new(polys, points, matrix, 3, initial_claims, 1, 1);
}


void test_lincheck_builder_new_with_valid_inputs() {
    const size_t N = 2, M = 2, NUM_VARS = 3, NUM_ACTIVE_VARS = 2;

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, 1 << NUM_VARS, f128_from_uint64(1));
    Points *points = points_init(NUM_VARS, f128_from_uint64(2) );
    Points* matrix_points = points_init(N * M * (1 << (NUM_ACTIVE_VARS * 2)), f128_from_uint64(1));
    MatrixLinear* matrix = matrix_linear_new(
        N * (1 << NUM_ACTIVE_VARS),
        M * (1 << NUM_ACTIVE_VARS),
        matrix_points->elems,
        matrix_points->len
    );
    F128 initial_claims[M] = { f128_from_uint64(4), f128_from_uint64(5) };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, NUM_ACTIVE_VARS, initial_claims, N, M);

    assert(builder->num_vars == NUM_VARS);
    assert(builder->num_active_vars == NUM_ACTIVE_VARS);

    // assert(builder->matrix == matrix);
    // assert(builder->points == points);
    //make sure the contents of builder->matrix are equal to matrix
    TEST_ASSERT_TRUE_MESSAGE(
        matrix_linear_equal(builder->matrix, matrix),
        "MatrixLinear objects are not equal"
    );
    TEST_ASSERT_TRUE_MESSAGE(
        points_equal(builder->points, points),
        "Points objects are not equal"
    );

    TEST_ASSERT_TRUE_MESSAGE(f128_eq(builder->initial_claims[0], f128_from_uint64(4)), 
        "Initial claim 0 mismatch");
    TEST_ASSERT_TRUE_MESSAGE(f128_eq(builder->initial_claims[1], f128_from_uint64(5)), 
        "Initial claim 1 mismatch");

    lincheck_builder_free(builder);
    mle_sequence_free(polys);
}

void test_lincheck_builder_build() {
    const size_t N = 2, M = 2, NUM_VARS = 3, NUM_ACTIVE_VARS = 2;

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, 1 << NUM_VARS, f128_zero());
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < (1 << NUM_VARS); j++) {
            polys->mle_poly[i].coeffs[j] = f128_from_uint64(i == 0 ? 5678 : 910);
        }
    }

    Points* points = points_init(NUM_VARS, f128_from_uint64(2));
    Points* matrix_points = points_init(N * M * (1 << (NUM_ACTIVE_VARS * 2)), f128_from_uint64(18));
    MatrixLinear *matrix = matrix_linear_new(
        N * (1 << NUM_ACTIVE_VARS),
        M * (1 << NUM_ACTIVE_VARS),
        matrix_points->elems, matrix_points->len
    );

    F128 initial_claims[M] = { f128_from_uint64(4), f128_from_uint64(5) };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, NUM_ACTIVE_VARS, initial_claims, N, M);
    F128 gamma = f128_from_uint64(1234);
    ProdCheck* prodcheck = lincheck_builder_build(builder, &gamma);

    TEST_ASSERT(prodcheck != NULL);
    // Optional: verify prodcheck->claim or contents of p_polys and q_polys

    prodcheck_free(prodcheck);
    lincheck_builder_free(builder);
    mle_sequence_free(polys);
}
