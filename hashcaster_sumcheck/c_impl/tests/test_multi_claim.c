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


void test_new_default_inputs(void) {
    // Polynomial: [1, 2]
    MLE_POLY poly = {
        .len = 2,
        .coeffs = malloc(2 * sizeof(F128))
    };
    poly.coeffs[0] = f128_from_uint64(1);
    poly.coeffs[1] = f128_from_uint64(2);

    // Points: [0]
    Points* points = points_init(1, f128_zero());

    // Openings: [0; 128]
    F128* openings = calloc(128, sizeof(F128));

    // Gamma powers: [1, 0, ..., 0]
    Points* gamma_pows = points_init(129, f128_zero());
    gamma_pows->elems[0] = f128_from_uint64(1);
    // F128* gamma_pows = calloc(129, sizeof(F128));
    // gamma_pows[0] = f128_from_uint64(1);

    // Polynomials: 2 identical MLEs: [1, 2, 3]
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 3, f128_zero());
    for (int i = 0; i < 2; ++i) {
        polys->mle_poly[i].coeffs[0] = f128_from_uint64(1);
        polys->mle_poly[i].coeffs[1] = f128_from_uint64(2);
        polys->mle_poly[i].coeffs[2] = f128_from_uint64(3);
    }

    MultiClaim* claim = multi_claim_new(&poly, points, openings, gamma_pows, polys);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            TEST_ASSERT_TRUE(f128_eq(claim->polys->mle_poly[i].coeffs[j], polys->mle_poly[i].coeffs[j]));
        }
    }

    TEST_ASSERT_TRUE(f128_eq(claim->gamma, f128_zero()));
    TEST_ASSERT_TRUE(f128_eq(claim->object->claim, f128_zero()));

    // Cleanup if necessary
}

void test_new_with_nonzero_openings(void) {
    // Create a non-default MultilinearLagrangianPolynomial
    MLE_POLY poly;
    poly.len = 2;
    poly.coeffs = malloc(2 * sizeof(F128));
    poly.coeffs[0] = f128_from_uint64(1);
    poly.coeffs[1] = f128_from_uint64(3);

    // Create points for evaluation
    Points* points = points_init(1, f128_zero());

    // Create openings with some non-zero values (all elements = 5)
    F128* openings = malloc(128 * sizeof(F128));
    for (int i = 0; i < 128; i++) {
        openings[i] = f128_from_uint64(5);
    }

    // Create gamma powers (0, 1, 2, ..., 128)
    Points* gamma_pows = points_init(129, f128_zero());
    for (int i = 0; i < 129; ++i) {
        gamma_pows->elems[i] = f128_from_uint64(i);
    }

    // Create polynomials (2 copies of [2, 4, 6])
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 3, f128_zero());
    for (int i = 0; i < 2; ++i) {
        polys->mle_poly[i].coeffs[0] = f128_from_uint64(2);
        polys->mle_poly[i].coeffs[1] = f128_from_uint64(4);
        polys->mle_poly[i].coeffs[2] = f128_from_uint64(6);
    }

    // Call the `new` function
    MultiClaim* claim = multi_claim_new(&poly, points, openings, gamma_pows, polys);

    // Compute expected initial claim
    F128 expected_claim = f128_zero();
    for (int i = 0; i < 128; i++) {
        F128 term = f128_mul(f128_from_uint64(i), f128_from_uint64(5));
        expected_claim = f128_add(expected_claim, term);
    }

    // Validate that claim->polys match the input polys
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            TEST_ASSERT_TRUE(f128_eq(claim->polys->mle_poly[i].coeffs[j], polys->mle_poly[i].coeffs[j]));
        }
    }

    // Validate gamma = 128
    TEST_ASSERT_TRUE(f128_eq(claim->gamma, f128_from_uint64(128)));

    // Validate the computed claim
    TEST_ASSERT_TRUE(f128_eq(claim->object->claim, expected_claim));

    // (Optional) free memory if needed
}


void test_multiclaim_complete(void) {
    const size_t NUM_VARS = 20;

    // Create a random multilinear polynomial with 2^NUM_VARS coefficients
    MLE_POLY* poly = mle_poly_random(1 << NUM_VARS);

    // Create random evaluation points
    Points* points = points_init(NUM_VARS, f128_zero());
    for (size_t i = 0; i < NUM_VARS; ++i) {
        points->elems[i] = f128_rand();
    }

    // Map the points to the inverse Frobenius orbit
    INVERSE_ORBIT_POINTS* points_inv_orbit = to_f128_inv_orbit(points);
    // Points* points_inv_orbit[128];
    // for (int i = 0; i < 128; ++i) {
    //     points_inv_orbit[i] = to_f128_inv_orbit
    // }

    // Evaluate the polynomial at the points in the orbit
    Evaluations* evaluations_inv_orbit = points_init(128, f128_zero());
    for (int i = 0; i < 128; ++i) {
        evaluations_inv_orbit->elems[i] =  mle_poly_evaluate_at(poly, &(points_inv_orbit->array_of_points[i]));
    }

    // Setup multiclaim builder
    // MultilinearLagrangianPolynomial polys[1];
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(1, poly->len, f128_zero());
    polys->mle_poly[0].len = poly->len;
    polys->mle_poly[0].coeffs = malloc(poly->len * sizeof(F128));
    for (size_t i = 0; i < poly->len; ++i) {
        polys->mle_poly[0].coeffs[i] = poly->coeffs[i];
    }
    MulticlaimBuilder* prover_builder = multiclaim_builder_new(polys, points, evaluations_inv_orbit);

    // Generate random gamma
    F128 gamma = f128_rand();

    // Build the prover
    MultiClaim* prover = multiclaim_builder_build(prover_builder, gamma);

    // Compute gamma powers
    Points* gamma_pows = compute_gammas_folding(gamma, 128);

    // Compute initial claim
    F128 claim = f128_zero();
    for (int i = 0; i < 128; ++i) {
        claim = f128_add (claim, f128_mul(gamma_pows->elems[i], evaluations_inv_orbit->elems[i]));
    }

    // Setup empty challenges
    Points* challenges = points_init(NUM_VARS, f128_zero());

    // Prover main loop
    for (size_t round = 0; round < NUM_VARS; ++round) {
        // Get round polynomial
        CompressedPoly* round_poly = prodcheck_round_polynomial(prover->object);
        UnivariatePolynomial* uncompressed_round_poly = uncompress_poly(round_poly);

        // Check that it has 3 coefficients
        TEST_ASSERT_EQUAL_UINT(3, uncompressed_round_poly->len);

        // Generate random challenge
        F128 challenge = f128_rand();

        // Update claim
        claim = univariate_polynomial_evaluate_at(uncompressed_round_poly, challenge);

        // // Add challenge to list
        // points_push(challenges, challenge);
        challenges->elems[round] = challenge;

        // Bind prover to challenge
        multi_claim_bind(prover, challenge);
    }

    // Compute equality evaluations
    F128 eq_evaluations = f128_zero();
    for (int i = 0; i < 128; ++i) {
        F128 eval = points_eq_eval(&(points_inv_orbit->array_of_points[i]), challenges);
        eq_evaluations = f128_add( f128_mul(gamma_pows->elems[i], eval), eq_evaluations);
    }

    // Compute expected final claim
    F128 expected = f128_mul(mle_poly_evaluate_at(poly, challenges), eq_evaluations);

    // Assert the final claim is correct
    TEST_ASSERT_TRUE(f128_eq(expected, claim));

    // Free memory (if necessary, depends on your memory model)
}