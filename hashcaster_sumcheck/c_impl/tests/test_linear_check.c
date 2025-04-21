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

// Assumes N is known and the same for both
bool prodcheck_equals(const ProdCheck* a, const ProdCheck* b) {
    if (a == NULL || b == NULL) return false;
    if (a->num_vars != b->num_vars) return false;
    if (!f128_eq(a->claim, b->claim)) return false;

    // Compare challenges
    if (!points_equal(a->challenges, b->challenges)) return false;

    // Compare polynomials
    for (size_t i = 0; i < a->N; i++) {
        if (!mle_poly_eq(&a->p_polys[i], &b->p_polys[i])) return false;
        if (!mle_poly_eq(&a->q_polys[i], &b->q_polys[i])) return false;
    }

    // Compare cached round msg (if both are NULL, they are equal)
    if ((a->cached_round_msg == NULL) != (b->cached_round_msg == NULL)) return false;

    if (a->cached_round_msg && b->cached_round_msg) {
        if (!compressed_poly_eq(a->cached_round_msg, b->cached_round_msg)) return false;
    }

    return true;
}

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
    // Setup constant values
    size_t chunk_size = 1 << NUM_ACTIVE_VARS;
    F128 val5678 = f128_from_raw(5678, 0);
    F128 val910 = f128_from_raw(910, 0);
    // F128 val2 = f128_from_raw(2, 0);
    // F128 val18 = f128_from_raw(18, 0);
    // F128 val4 = f128_from_raw(4, 0);
    // F128 val5 = f128_from_raw(5, 0);
    F128 gamma = f128_from_raw(1234, 0);
    ProdCheck* prodcheck = lincheck_builder_build(builder, &gamma);

    TEST_ASSERT(prodcheck != NULL);
    // Expected p_polys
    MLE_POLY *expected_p0 = mle_poly_from_constant(chunk_size,val5678);
    MLE_POLY *expected_p1 = mle_poly_from_constant(chunk_size,val910);

    // 18692268690725379462709786785192376188 = high: 0x0e10000000000000, low: 0x000000000000777c
    F128 qval = f128_from_raw(0x000000000000777c, 0x0e10000000000000);
    MLE_POLY *expected_q0 = mle_poly_from_constant(chunk_size,qval);
    MLE_POLY *expected_q1 = mle_poly_from_constant(chunk_size,qval);

    // 258410230055561301416705741312625744282 = high: 0x0399000000000000, low: 0x000000000000abcd
    F128 expected_claim = f128_from_raw(0x000000000000199a, 0xc268000000000000);

    // expected mle_polys
    MLE_POLY expected_p_arr[2] = { *expected_p0, *expected_p1 };
    MLE_POLY expected_q_arr[2] = { *expected_q0, *expected_q1 };

    // Expected ProdCheck
    ProdCheck *expected = prodcheck_new(expected_p_arr, expected_q_arr, N, expected_claim, 0);

    // Compare
    TEST_ASSERT_TRUE(prodcheck_equals(prodcheck, expected));    // Optional: verify prodcheck->claim or contents of p_polys and q_polys

    prodcheck_free(prodcheck);
    prodcheck_free(expected);
    // Free expected polynomials
    mle_poly_free(expected_p0);
    mle_poly_free(expected_p1);
    mle_poly_free(expected_q0);
    mle_poly_free(expected_q1);
    // Free builder resources
    // matrix_linear_free(builder->matrix);
    // points_free(builder->points);
    // points_free(matrix_points);
    // points_free(points);
    lincheck_builder_free(builder);
    mle_sequence_free(polys);
}


void test_lincheck_builder_invalid_matrix_input_size(void) {
    const size_t N = 2, M = 2, NUM_VARS = 3, NUM_ACTIVE_VARS = 2;

    size_t n_in = (N - 1) * (1 << NUM_ACTIVE_VARS); // incorrect
    size_t n_out = M * (1 << NUM_ACTIVE_VARS);
    size_t len = n_in * n_out;

    Points* matrix_points = points_init(len, f128_from_uint64(1));
    MatrixLinear* matrix = matrix_linear_new(n_in, n_out, matrix_points->elems, matrix_points->len);

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, 1 << NUM_VARS, f128_from_uint64(1));
    Points* points = points_init(NUM_VARS, f128_from_uint64(2));
    F128 initial_claims[M] = { f128_from_uint64(4), f128_from_uint64(5) };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, NUM_ACTIVE_VARS, initial_claims, N, M);

    TEST_ASSERT_NULL(builder); // should fail due to invalid matrix input size

    mle_sequence_free(polys);
    matrix_linear_free(matrix);
    points_free(points);
    points_free(matrix_points);
}

void test_lincheck_builder_invalid_matrix_output_size(void) {
    const size_t N = 2, M = 2, NUM_VARS = 3, NUM_ACTIVE_VARS = 2;

    size_t n_in = N * (1 << NUM_ACTIVE_VARS);
    size_t n_out = (M + 1) * (1 << NUM_ACTIVE_VARS); // invalid output size
    size_t len = n_in * n_out;

    Points* matrix_points = points_init(len, f128_from_uint64(1));
    MatrixLinear* matrix = matrix_linear_new(n_in, n_out, matrix_points->elems, matrix_points->len);

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, 1 << NUM_VARS, f128_from_uint64(1));
    Points* points = points_init(NUM_VARS, f128_from_uint64(2));
    F128 initial_claims[M] = { f128_from_uint64(4), f128_from_uint64(5) };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, NUM_ACTIVE_VARS, initial_claims, N, M);

    TEST_ASSERT_NULL(builder); // should fail

    mle_sequence_free(polys);
    matrix_linear_free(matrix);
    points_free(points);
    points_free(matrix_points);
}

void test_lincheck_builder_insufficient_points(void) {
    const size_t N = 2, M = 2, NUM_VARS = 2, NUM_ACTIVE_VARS = 3;

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, 1 << NUM_VARS, f128_from_uint64(1));
    Points* points = points_init(NUM_VARS, f128_from_uint64(2));
    Points* matrix_points = points_init(N * M * (1 << (NUM_ACTIVE_VARS * 2)), f128_from_uint64(1));
    MatrixLinear* matrix = matrix_linear_new(
        N * (1 << NUM_ACTIVE_VARS),
        M * (1 << NUM_ACTIVE_VARS),
        matrix_points->elems, matrix_points->len
    );

    F128 initial_claims[M] = { f128_from_uint64(4), f128_from_uint64(5) };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, NUM_ACTIVE_VARS, initial_claims, N, M);

    TEST_ASSERT_NULL(builder); // not enough vars to satisfy active ones

    mle_sequence_free(polys);
    matrix_linear_free(matrix);
    points_free(points);
    points_free(matrix_points);
}

void test_lincheck_builder_invalid_polynomial_length(void) {
    const size_t N = 2, M = 2, NUM_VARS = 3, NUM_ACTIVE_VARS = 2;

    // Off-by-one polynomial size
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, (1 << NUM_VARS) - 1, f128_from_uint64(1));
    Points* points = points_init(NUM_VARS, f128_from_uint64(2));
    Points* matrix_points = points_init(N * M * (1 << (NUM_ACTIVE_VARS * 2)), f128_from_uint64(1));
    MatrixLinear* matrix = matrix_linear_new(
        N * (1 << NUM_ACTIVE_VARS),
        M * (1 << NUM_ACTIVE_VARS),
        matrix_points->elems, matrix_points->len
    );

    F128 initial_claims[M] = { f128_from_uint64(4), f128_from_uint64(5) };

    LinCheckBuilder* builder = lincheck_builder_new(polys, points, matrix, NUM_ACTIVE_VARS, initial_claims, N, M);

    TEST_ASSERT_NULL(builder); // invalid poly length

    mle_sequence_free(polys);
    matrix_linear_free(matrix);
    points_free(points);
    points_free(matrix_points);
}