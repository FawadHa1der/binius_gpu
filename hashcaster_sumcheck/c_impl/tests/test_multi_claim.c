#include "test_multi_claim.h"
#include "unity.h"
#include "../multi_claim.h"
#include "../prod_check.h"
#include "../field.h"
#include "../mle_poly.h"
#include "../univariate_poly.h"
#include "../linear_check.h"
#include "../efficient_matrix.h"
#include "../evaluations.h"

// void test_multiclaim_builder_default(void) {
//     MLE_POLY_SEQUENCE* polys = mle_sequence_new(3, 0, f128_zero());
//     // for (int i = 0; i < 3; ++i) {
//     //     polys[i] = mle_
//     // }

//     Points points = points_default();
//     Evaluations* openings = points_init(3 * 128, f128_zero());

//     MulticlaimBuilder* builder = multiclaim_builder_new(polys, &points, openings);

//     for (int i = 0; i < 3; ++i) {
//         TEST_ASSERT_TRUE(builder->polys->mle_poly[i].len == 0);
//     }

//     TEST_ASSERT_TRUE(builder->points->len == 0);
//     TEST_ASSERT_TRUE(openings->len == 0);
// }


void test_multiclaim_builder_new_valid(void) {
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 2, f128_zero());
    polys->mle_poly->coeffs[0] = f128_from_uint64(1);
    polys->mle_poly->coeffs[1] = f128_from_uint64(2);
    polys->mle_poly[1].coeffs[0] = f128_from_uint64(3);
    polys->mle_poly[1].coeffs[1] = f128_from_uint64(4);

    Points* points = points_init(1,f128_from_uint64(1) );
    Evaluations* openings = points_init(2 * 128, f128_zero());

    MulticlaimBuilder* builder = multiclaim_builder_new(polys, points, openings);

    TEST_ASSERT_EQUAL_UINT(builder->polys[0].len, 2);
    TEST_ASSERT_TRUE(points_equal(builder->points, points));
    TEST_ASSERT_TRUE(points_equal(builder->openings, openings));
}

void test_multiclaim_builder_new_invalid_polynomial_length(void) {
    MLE_POLY_SEQUENCE* polys = malloc(sizeof(MLE_POLY_SEQUENCE));
    polys->len = 2;
    polys->mle_poly = malloc(2 * sizeof(MLE_POLY));
    polys->mle_poly[0].len = 1;
    polys->mle_poly[0].coeffs = malloc(sizeof(F128));
    polys->mle_poly[0].coeffs[0] = f128_from_uint64(1);
    polys->mle_poly[1].len = 2;
    polys->mle_poly[1].coeffs = malloc(2 * sizeof(F128));
    polys->mle_poly[1].coeffs[0] = f128_from_uint64(3);
    polys->mle_poly[1].coeffs[1] = f128_from_uint64(4);

    Points* points = points_init(1, f128_from_uint64(1));
    Evaluations* openings = points_init(2 * 128, f128_zero());
    multiclaim_builder_new(polys, points, openings);

    bool failed = false;
    TEST_ASSERT_TRUE_MESSAGE(failed, "Expected validation failure due to incorrect polynomial length");
}

void test_multiclaim_builder_new_edge_case(void) {
    MLE_POLY_SEQUENCE* polys = malloc(sizeof(MLE_POLY_SEQUENCE));
    polys->len = 1;
    polys->mle_poly = mle_poly_from_constant(1, f128_from_uint64(42));

    // MLE_POLY* polys[1] = {
    //     mlag_poly_from_vals((F128[]){f128_from_uint64(42)}, 1)
    // };
    Points* points = points_init(0, f128_zero());
    Evaluations* openings = points_init(128, f128_zero());

    MulticlaimBuilder* builder = multiclaim_builder_new(polys, points, openings);

    TEST_ASSERT_EQUAL_UINT(builder->polys->mle_poly[0].len, 1);
    TEST_ASSERT_TRUE(builder->points->len == points->len);
    TEST_ASSERT_TRUE(points_equal(builder->openings, openings));
}


void test_multiclaim_builder_build_simple_case(void) {
    F128 p11 = f128_from_uint64(1), p12 = f128_from_uint64(2);
    F128 p21 = f128_from_uint64(3), p22 = f128_from_uint64(4);

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 2, f128_from_uint64(2));
    polys->mle_poly[0].coeffs[0] = p11;
    polys->mle_poly[0].coeffs[1] = p12;
    polys->mle_poly[1].coeffs[0] = p21;
    polys->mle_poly[1].coeffs[1] = p22;

    Points* points = points_init(1, f128_from_uint64(1));
    Evaluations* openings = points_init(2 * 128, f128_zero());
    for (int i = 0; i < 2 * 128; ++i) {
        openings->elems[i] = f128_from_uint64(i);
    }

    MulticlaimBuilder* builder = multiclaim_builder_new(polys, points, openings);
    F128 gamma = f128_from_uint64(2);

    MultiClaim* claim = multiclaim_builder_build(builder, gamma);

    Points *gamma_pows = compute_gammas_folding(gamma, 256);

    F128 expected_poly[2] = {
        f128_add(p11, f128_mul(p21, gamma_pows->elems[128])),
        f128_add(p12, f128_mul(p22, gamma_pows->elems[128]))
    };

    TEST_ASSERT_TRUE(f128_eq(claim->object->p_polys[0].coeffs[0], expected_poly[0]));
    TEST_ASSERT_TRUE(f128_eq(claim->object->p_polys[0].coeffs[1], expected_poly[1]));

    F128 expected_claim = f128_zero();
    for (int i = 0; i < 128; ++i) {
        F128 temp = f128_add(openings->elems[i], f128_mul(openings->elems[i + 128], gamma_pows->elems[128]));
        expected_claim = f128_add(expected_claim, f128_mul(temp, gamma_pows->elems[i]));
    }

    TEST_ASSERT_TRUE(f128_eq(claim->object->claim, expected_claim));
}