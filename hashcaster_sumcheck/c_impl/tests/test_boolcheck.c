#include "test_boolcheck.h"
#include "unity.h"
#include "time.h"
#include "../multi_claim.h"


void test_trit_mapping_small_c(void) {
    uint16_t expected_bit[] = {0, 1, 3, 4};
    uint16_t expected_trit[] = {0, 2, 1, 4, 6, 1, 3, 3, 3};
    size_t bit_len = 0, trit_len = 0;
    uint16_t *bit_mapping = NULL;
    uint16_t *trit_mapping = NULL;

    compute_trit_mappings(1, &bit_mapping, &bit_len, &trit_mapping, &trit_len);

    TEST_ASSERT_EQUAL_UINT32(4, bit_len);
    TEST_ASSERT_EQUAL_UINT32(9, trit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_bit, bit_mapping, bit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_trit, trit_mapping, trit_len);

    free(bit_mapping);
    free(trit_mapping);
}

void test_trit_mapping_medium_c(void) {
    uint16_t expected_bit[] = {0, 1, 3, 4, 9, 10, 12, 13};
    uint16_t expected_trit[] = {
        0, 2, 1, 4, 6, 1, 3, 3, 3, 8, 10, 1, 12, 14, 1, 3, 3, 3,
        9, 9, 9, 9, 9, 9, 9, 9, 9
    };
    size_t bit_len = 0, trit_len = 0;
    uint16_t *bit_mapping = NULL;
    uint16_t *trit_mapping = NULL;

    compute_trit_mappings(2, &bit_mapping, &bit_len, &trit_mapping, &trit_len);

    TEST_ASSERT_EQUAL_UINT32(8, bit_len);
    TEST_ASSERT_EQUAL_UINT32(27, trit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_bit, bit_mapping, bit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_trit, trit_mapping, trit_len);

    free(bit_mapping);
    free(trit_mapping);
}

void test_trit_mapping_large_c(void) {
    uint16_t expected_bit[] = {0, 1, 3, 4, 9, 10, 12, 13, 27, 28, 30, 31, 36, 37, 39, 40, 81, 82, 84, 85, 90, 91,
        93, 94, 108, 109, 111, 112, 117, 118, 120, 121
};
    uint16_t expected_trit[] = {
        0, 2, 1, 4, 6, 1, 3, 3, 3, 8, 10, 1, 12, 14, 1, 3, 3, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        16, 18, 1, 20, 22, 1, 3, 3, 3, 24, 26, 1, 28, 30, 1, 3, 3, 3, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
        27, 27, 27, 27, 27, 27, 27, 27, 32, 34, 1, 36, 38, 1, 3, 3, 3, 40, 42, 1, 44, 46,
        1, 3, 3, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 48, 50, 1, 52, 54, 1, 3, 3, 3, 56, 58, 1,
        60, 62, 1, 3, 3, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 27, 27, 27, 27, 27, 27, 27, 27, 27,
        27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 81, 81, 81,
        81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81,
        81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81,
        81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81,
        81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81
    };
    size_t bit_len = 0, trit_len = 0;
    uint16_t *bit_mapping = NULL;
    uint16_t *trit_mapping = NULL;

    compute_trit_mappings(4, &bit_mapping, &bit_len, &trit_mapping, &trit_len);

    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_bit, bit_mapping, bit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_trit, trit_mapping, trit_len);

    free(bit_mapping);
    free(trit_mapping);
}



void test_trit_mapping_no_c(void) {
    uint16_t expected_bit[] = {0, 1};
    uint16_t expected_trit[] = {0, 2, 1};
    size_t bit_len = 0, trit_len = 0;
    uint16_t *bit_mapping = NULL;
    uint16_t *trit_mapping = NULL;

    compute_trit_mappings(0, &bit_mapping, &bit_len, &trit_mapping, &trit_len);

    TEST_ASSERT_EQUAL_UINT32(2, bit_len);
    TEST_ASSERT_EQUAL_UINT32(3, trit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_bit, bit_mapping, bit_len);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected_trit, trit_mapping, trit_len);

    free(bit_mapping);
    free(trit_mapping);
}

void test_package_algebraic(
    const Algebraic_Params* params,
    const Points* data,
    size_t idx_a,
    size_t offset,
    F128 *ret[3] // TODO a better data structure for a 2D array ?
) {

    // not being used by the extend_n_tables function
}

void test_package_linear(
    const Algebraic_Params* params,
    const Points* data,
    Points* out
) {
    F128 res = f128_zero();
    for (size_t i = 0; i < params->input_size; i++) {
        res = f128_add(res,data->elems[i] ) ;
    }

    for (size_t i = 0; i < params->output_size; i++) {
        out->elems[i] = res;
    }
}

void test_package_quadratic(
    const Algebraic_Params* params,
    const Points* data,
    Points* out
) {
    F128 res = f128_zero();
    for (size_t i = 0; i < params->input_size; i++) {
        res = f128_add(res, f128_mul(data->elems[i], data->elems[i]) );
    }
    for (size_t i = 0; i < params->output_size; i++) {
        out->elems[i] = res;
    }
}

void test_extend_n_tables(void) {
    const size_t dims = 3;
    const size_t N = 3;
    const size_t M = 1;
    const size_t C = 2;
    const size_t table_size = 1 << dims;

    // Initialize tables
    // F128 table1[table_size], table2[table_size], table3[table_size];
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(N, table_size, f128_zero());
    for (size_t i = 0; i < N; ++i) {
        polys->mle_poly[i].len = table_size;
        polys->mle_poly[i].coeffs = (F128*)malloc(sizeof(F128) * table_size);
        for (size_t j = 0; j < table_size; ++j) {
            polys->mle_poly[i].coeffs[j] = f128_from_uint64(j + (i* 10));
        }
    }
   
    for (size_t i = 0; i < table_size; ++i) {
        // table1[i] = f128_from_uint64(i);
        // table2[i] = f128_from_uint64(i + 10);
        // table3[i] = f128_from_uint64(i + 20);
    }
    // const F128* tables[N] = {table1, table2, table3};

    // Dummy algebraic logic
    // DummyAlgebraicPackage dummy = get_dummy_package();
    Algebraic_Params dummy_params;
    dummy_params.input_size = 3;
    dummy_params.output_size = 1;
    Algebraic_Functions dummy_funcs;
    dummy_funcs.algebraic = test_package_algebraic;
    dummy_funcs.linear = test_package_linear;
    dummy_funcs.quadratic = test_package_quadratic;


    // Setup points and claims
    Points* pts = points_init(dims, f128_zero());
    Points* claims = points_init(M, f128_zero());
    Points* gammas = points_init(M, f128_one());

    // Setup polynomial sequence
    // MLE_POLY_SEQUENCE* polys = mle_sequence_from_tables((F128**)tables, N, table_size);

    // BoolCheckBuilder
    BoolCheckBuilder* builder = bool_check_builder_new(C,
        pts, claims, polys, &dummy_params, &dummy_funcs);
    builder ->gammas = gammas;        

    // Compute trit mapping
    uint16_t *bit_map = NULL, *trit_map = NULL;
    size_t bit_len = 0, trit_len = 0;
    compute_trit_mappings(C, &bit_map, &bit_len, &trit_map, &trit_len);

    // Extend tables
    Points* result = extend_n_tables(polys, N, dims, C, trit_map, builder,
                                   linear_compressed, quadratic_compressed);

    // Expected output
    F128 expected[] = {
        {0x00000000000001b7ULL, 0xd750000000000000ULL},
        {0x00000000000001b7ULL, 0x4554000000000000ULL},
        {0x0000000000000001ULL, 0x9204000000000000ULL},
        {0x00000000000001aaULL, 0xe100000000000000ULL},
        {0x00000000000001aaULL, 0x7304000000000000ULL},
        {0x0000000000000001ULL, 0x9204000000000000ULL},
        {0x000000000000001bULL, 0x3650000000000000ULL},
        {0x000000000000001bULL, 0x3650000000000000ULL},
        {0x0000000000000000ULL, 0x0000000000000000ULL},
        {0x00000000000001d7ULL, 0x0e10000000000000ULL},
        {0x00000000000001d7ULL, 0x9c14000000000000ULL},
        {0x0000000000000001ULL, 0x9204000000000000ULL},
        {0x0000000000000060ULL, 0xd940000000000000ULL},
        {0x0000000000000060ULL, 0x4b44000000000000ULL},
        {0x0000000000000001ULL, 0x9204000000000000ULL},
        {0x00000000000001a9ULL, 0xd750000000000000ULL},
        {0x00000000000001a9ULL, 0xd750000000000000ULL},
        {0x0000000000000000ULL, 0x0000000000000000ULL},
        {0x000000000000006cULL, 0xd940000000000000ULL},
        {0x000000000000006cULL, 0xd940000000000000ULL},
        {0x0000000000000000ULL, 0x0000000000000000ULL},
        {0x00000000000001deULL, 0x3840000000000000ULL},
        {0x00000000000001deULL, 0x3840000000000000ULL},
        {0x0000000000000000ULL, 0x0000000000000000ULL},
        {0x00000000000001b2ULL, 0xe100000000000000ULL},
        {0x00000000000001b2ULL, 0xe100000000000000ULL},
        {0x0000000000000000ULL, 0x0000000000000000ULL},
    };

    for (size_t i = 0; i < 27; ++i) {
        TEST_ASSERT_TRUE(f128_eq(result->elems[i], expected[i]));
    }

    // Cleanup
    free(bit_map);
    free(trit_map);
    free(result);
    points_free(pts);
    points_free(claims);
    mle_sequence_free(polys);
    bool_check_builder_free(builder);
}

void test_new_andcheck1(void) {
    const size_t PHASE_SWITCH = 5;
    const size_t num_vars = 20;

    Points * points = points_init(num_vars, f128_zero());
    for (size_t i = 0; i < num_vars; ++i) {
        points->elems[i] = f128_rand();
    }

    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 1 << num_vars, f128_zero());
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < (1 << num_vars); ++j) {
            polys->mle_poly[i].coeffs[j] = f128_rand();
        }
    }

    MLE_POLY* p = &polys->mle_poly[0];
    MLE_POLY* q = &polys->mle_poly[1];
    MLE_POLY* pq = mle_poly_from_constant(1 << num_vars, f128_zero());
    for (size_t i = 0; i < (1 << num_vars); ++i) {
        pq->coeffs[i] = f128_bitand(p->coeffs[i], q->coeffs[i] );
    }

    F128 initial_claim = mle_poly_evaluate_at(pq, points);
    F128 current_claim = initial_claim;
    F128 gamma = f128_rand();
    Points* challenges = points_init(0, f128_zero());

    Algebraic_Params dummy_params;
    dummy_params.input_size = 2;
    dummy_params.output_size = 1;

    Algebraic_Functions dummy_funcs;
    dummy_funcs.algebraic = and_package_algebraic;
    dummy_funcs.linear = and_package_linear;
    dummy_funcs.quadratic = and_package_quadratic;

    BoolCheckBuilder *builder = bool_check_builder_new(
        PHASE_SWITCH, points, points_init(1, initial_claim) , polys, &dummy_params, &dummy_funcs);
    BoolCheck *boolcheck = boolcheck_new(builder, gamma);
    
    for (size_t i = 0; i < num_vars; ++i) {
        CompressedPoly *round_poly = boolcheck_round_polynomial(boolcheck);
        F128 r = f128_rand();
        UnivariatePolynomial* coeffs = uncompress_poly(round_poly, current_claim );
        // current_claim = fixed_univar_evaluate(&coeffs, r);
        current_claim = univariate_polynomial_evaluate_at(coeffs, r);
        boolcheck_bind(boolcheck, &r);
        points_push(challenges, r);
    }

    BoolCheckOutput *out = boolcheck_finish(boolcheck);
    untwist_evals(out->frob_evals);

    F128* and_algebraic_output[3];
    for (int i = 0; i < 3; i++) {
        and_algebraic_output[i] = (F128*)malloc(sizeof(F128) * 1);
        and_algebraic_output[i][0] = f128_zero();
    }
    and_package_algebraic(
        &dummy_params, out->frob_evals, 0, 1, and_algebraic_output
    );

    F128 expected = and_algebraic_output[0][0];
    expected = f128_mul(expected, points_eq_eval(points, challenges));

    TEST_ASSERT_TRUE(f128_eq(current_claim, expected));
}

void test_new_andcheck_with_multiclaim(void) {
    const size_t PHASE_SWITCH = 5;
    const size_t num_vars = 20;

    // Setup random points
    Points *points = points_init(num_vars, f128_zero());
    for (size_t i = 0; i < num_vars; ++i) {
        points->elems[i] = f128_from_uint64(i);
    }

    // Generate two random multilinear polynomials p and q
    MLE_POLY_SEQUENCE* polys = mle_sequence_new(2, 1 << num_vars, f128_zero());
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < (1 << num_vars); ++j) {
            polys->mle_poly[i].coeffs[j] = f128_from_uint64(1);
        }
    }
    MLE_POLY* p = &polys->mle_poly[0];
    MLE_POLY* q = &polys->mle_poly[1];

    // Compute their AND
    MLE_POLY* pq = mle_poly_from_constant(1 << num_vars, f128_zero());
    for (size_t i = 0; i < (1 << num_vars); ++i) {
        pq->coeffs[i] = f128_bitand(p->coeffs[i], q->coeffs[i]);
    }

    // Evaluate AND at points to get initial_claim
    F128 initial_claim = mle_poly_evaluate_at(pq, points);
    F128 current_claim = initial_claim;
    Points* challenges = points_init(0, f128_zero());

    // Setup algebraic package
    Algebraic_Params and_params;
    and_params.input_size = 2;
    and_params.output_size = 1;
    Algebraic_Functions and_funcs;
    and_funcs.algebraic = and_package_algebraic;
    and_funcs.linear = and_package_linear;
    and_funcs.quadratic = and_package_quadratic;

    // Setup random gamma
    F128 gamma = f128_from_uint64(10);

    // Setup BoolCheckBuilder and BoolCheck
    BoolCheckBuilder *builder = bool_check_builder_new(
        PHASE_SWITCH, points, points_init(1, initial_claim), polys, &and_params, &and_funcs);
    BoolCheck *boolcheck = boolcheck_new(builder, gamma);

    // Timing
    clock_t t_start = clock();

    // Main boolcheck rounds
    for (size_t i = 0; i < num_vars; ++i) {
        CompressedPoly *round_poly = boolcheck_round_polynomial(boolcheck);
        F128 r = f128_from_uint64(i);
        UnivariatePolynomial* coeffs = uncompress_poly(round_poly, current_claim);
        current_claim = univariate_polynomial_evaluate_at(coeffs, r);
        boolcheck_bind(boolcheck, &r);
        points_push(challenges, r);
        free(coeffs);
    }

    BoolCheckOutput *out = boolcheck_finish(boolcheck);
    Evaluations* frob_evals_copy = points_copy( out->frob_evals);

    // Untwist the output Frobenius evaluations
    untwist_evals(frob_evals_copy);

    // Compute expected value using and_package_algebraic
    F128* and_algebraic_output[3];
    for (int i = 0; i < 3; i++) {
        and_algebraic_output[i] = (F128*)malloc(sizeof(F128) * 1);
        and_algebraic_output[i][0] = f128_zero();
    }
    and_package_algebraic(&and_params, frob_evals_copy, 0, 1, and_algebraic_output);
    F128 expected = and_algebraic_output[0][0];
    expected = f128_mul(expected, points_eq_eval(points, challenges));

    TEST_ASSERT_TRUE(f128_eq(current_claim, expected));

    //////////////////////////////MULTI CALIM //////////////////////////////////////
    
    Points* challenges_copy = points_copy(challenges);
    INVERSE_ORBIT_POINTS* points_inv_orbit = to_f128_inv_orbit(challenges_copy);

    // Print timing
    printf("Execution time of Boolcheck: %ld ms\n", (clock() - t_start) * 1000 / CLOCKS_PER_SEC);

    // --- Multiclaim part ---
    t_start = clock();
    
    
    MulticlaimBuilder* prover_builder = multiclaim_builder_new(polys, challenges_copy, out->frob_evals);

    // // Setup multiclaim builder with p, q, and the Frobenius evaluations from BoolCheck
    // MulticlaimBuilder* multiclaim_builder = multiclaim_builder_new(
    // );

    // New gamma for folding
    F128 gamma2 = f128_from_uint64(20);
    F128 gamma128 = f128_pow(gamma2, 128);


    MultiClaim* multiclaim_prover = multiclaim_builder_build(prover_builder, gamma2);
    // Compute the claim as evaluation of frob_evals at gamma2
    Points* first_elem_frob_evals = points_init(1, out->frob_evals->elems[0]);
    F128 claim =  univariate_polynomial_evaluate_at(first_elem_frob_evals, gamma2);
    Points* multiclaim_challenges = points_init(0, f128_zero());

    for (size_t i = 0; i < num_vars; ++i) {
        CompressedPoly* round_poly =  prodcheck_round_polynomial(multiclaim_prover->object);
        UnivariatePolynomial* coeffs = uncompress_poly(round_poly, claim);
        F128 challenge = f128_from_uint64(i);
        claim = univariate_polynomial_evaluate_at(coeffs, challenge);
        points_push(multiclaim_challenges, challenge);

        multi_claim_bind(multiclaim_prover, challenge);
        free(coeffs);
    }

    Evaluations* multi_claim_finish_evaluations = multi_claim_finish(multiclaim_prover);

    // Finish multiclaim
    // Points* eq_evals_pts[128];
    Points* eq_evals_pts = points_init(128, f128_zero());
    for (size_t i = 0; i < 128; ++i) {
        // eq_evals_pts->elems[j] = multiclaim_challenges->elems[j];
        eq_evals_pts->elems[i] = points_eq_eval(&points_inv_orbit->array_of_points[i], multiclaim_challenges);
    }
    // Compute eq_evaluations at gamma2
    F128 eq_evaluation = univariate_polynomial_evaluate_at(eq_evals_pts, gamma2);

    // Validate the claim
    F128 output_eval = univariate_polynomial_evaluate_at(multi_claim_finish_evaluations, gamma128);
    F128 check = f128_mul(output_eval, eq_evaluation);
    TEST_ASSERT_TRUE(f128_eq(check, claim));

    // Check that p and q evaluated at challenges match the outputs
    F128 p_eval = mle_poly_evaluate_at(p, multiclaim_challenges);
    F128 q_eval = mle_poly_evaluate_at(q, multiclaim_challenges);

    TEST_ASSERT_TRUE(f128_eq(p_eval, multi_claim_finish_evaluations->elems[0]));
    TEST_ASSERT_TRUE(f128_eq(q_eval, multi_claim_finish_evaluations->elems[1]));

    printf("Execution time of Multiclaim: %ld ms\n", (clock() - t_start) * 1000 / CLOCKS_PER_SEC);

    // Cleanup
    for (int i = 0; i < 3; i++) free(and_algebraic_output[i]);
    points_free(points);
    points_free(challenges);
    points_free(challenges_copy);
    points_free(multiclaim_challenges);
    points_free(eq_evals_pts);
    mle_sequence_free(polys);
    mle_poly_free(pq);
    bool_check_builder_free(builder);
    // boolcheck_free(boolcheck);
    boolcheck_output_free(out);
    // multi_claim_builder_free(prover_builder);
    // multi_claim_free(multiclaim_prover);
}


void test_compute_imaginary_rounds(void) {
    // Set a phase switch parameter, which controls the folding phases.
    const size_t PHASE_SWITCH = 2;
    const size_t num_vars = 3;

    // Construct input points
    Points *points = points_init(num_vars, f128_zero());
    for (size_t i = 0; i < num_vars; ++i) {
        points->elems[i] = f128_from_uint64(i);
    }

    // Construct p and q as polynomials of length 1 << num_vars
    size_t poly_len = 1 << num_vars;
    MLE_POLY_SEQUENCE *polys = mle_sequence_new(2, poly_len, f128_zero());
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < poly_len; ++j) {
            polys->mle_poly[i].coeffs[j] = f128_from_uint64(j);
        }
    }

    // Compute AND of p and q
    MLE_POLY *p = &polys->mle_poly[0];
    MLE_POLY *q = &polys->mle_poly[1];
    MLE_POLY *pq = mle_poly_from_constant(poly_len, f128_zero());
    for (size_t i = 0; i < poly_len; ++i) {
        pq->coeffs[i] = f128_bitand(p->coeffs[i], q->coeffs[i]);
    }

    // Evaluate AND at points to get initial_claim
    F128 initial_claim = mle_poly_evaluate_at(pq, points);
    F128 current_claim = initial_claim;

    // Set gamma and r
    F128 gamma = f128_from_uint64(1234);
    F128 r = f128_from_uint64(5678);

    // Setup algebraic package
    Algebraic_Params and_params;
    and_params.input_size = 2;
    and_params.output_size = 1;
    Algebraic_Functions and_funcs;
    and_funcs.algebraic = and_package_algebraic;
    and_funcs.linear = and_package_linear;
    and_funcs.quadratic = and_package_quadratic;

    // Setup BoolCheckBuilder and BoolCheck
    BoolCheckBuilder *builder = bool_check_builder_new(
        PHASE_SWITCH, points, points_init(1, initial_claim), polys, &and_params, &and_funcs);
    BoolCheck *boolcheck = boolcheck_new(builder, gamma);

    // Validate initial extended_table
    F128 expected_table[] = {
        f128_zero(),
        f128_from_uint64(1),
        f128_from_uint64(1),
        f128_from_uint64(2),
        f128_from_uint64(3),
        f128_from_uint64(1),
        f128_from_uint64(2),
        f128_from_uint64(2),
        f128_zero(),
        f128_from_uint64(4),
        f128_from_uint64(5),
        f128_from_uint64(1),
        f128_from_uint64(6),
        f128_from_uint64(7),
        f128_from_uint64(1),
        f128_from_uint64(2),
        f128_from_uint64(2),
        f128_zero(),
        f128_from_uint64(4),
        f128_from_uint64(4),
        f128_zero(),
        f128_from_uint64(4),
        f128_from_uint64(4),
        f128_zero(),
        f128_zero(),
        f128_zero(),
        f128_zero()
    };
    for (size_t i = 0; i < 27; ++i) {
        TEST_ASSERT_TRUE(f128_eq(boolcheck->extended_table->elems[i], expected_table[i]));
    }

    // Compute the round polynomial for the imaginary round
    CompressedPoly *compressed_round_poly = boolcheck_round_polynomial(boolcheck);

    // Check compressed coefficients
    // These values should match the Rust test values
    F128 expected_compressed[3] = {
        (F128){.low = 0x000000000000000dULL, .high = 0xfa28000000000000}, // 332514690820570361331092984923254947853
        f128_from_uint64(1),
        f128_from_uint64(1)
    };
    for (int i = 0; i < 3; ++i) {
        TEST_ASSERT_TRUE(f128_eq(compressed_round_poly->compressed_coeff->elems[i], expected_compressed[i]));
    }

    // Decompress the polynomial and check coefficients
    UnivariatePolynomial *round_poly = uncompress_poly(compressed_round_poly, initial_claim);
    F128 expected_decompressed[4] = {
        (F128){.low = 0x000000000000000dULL, .high = 0xfa28000000000000},
        (F128){.low = 0x000000000000000dULL, .high = 0xfa28000000000000},
        f128_from_uint64(1),
        f128_from_uint64(1)
    };
    for (int i = 0; i < 4; ++i) {
        TEST_ASSERT_TRUE(f128_eq(round_poly->elems[i], expected_decompressed[i]));
    }

    // Evaluate the decompressed polynomial at r and verify the result
    F128 eval_result = univariate_polynomial_evaluate_at(round_poly, r);
    F128 expected_eval = (F128){.low = 0x0000001ebe1310a2, .high = 0xd5cb5d8000000000}; // 284181495769767592368287233794578256034
    TEST_ASSERT_TRUE(f128_eq(eval_result, expected_eval));

    // Bind r and verify updates
    boolcheck_bind(boolcheck, &r);

    // Validate updated extended_table (length 9)
    F128 expected_table2[] = {
        (F128){.low = 0x0000000001404154, .high = 0x6d0ea00000000000}, // 144961788882111902000582228079393390932
        (F128){.low = 0x0000000001404156, .high = 0x6d0ea00000000000}, // 144961788882111902000582228079393390934
        f128_from_uint64(2),
        (F128){.low = 0x0000000001404150, .high = 0x6d0ea00000000000}, // 144961788882111902000582228079393390928
        (F128){.low = 0x0000000001404152, .high = 0x6d0ea00000000000}, // 144961788882111902000582228079393390930
        f128_from_uint64(2),
        f128_from_uint64(4),
        f128_from_uint64(4),
        f128_zero()
    };
    for (size_t i = 0; i < 9; ++i) {
        TEST_ASSERT_TRUE(f128_eq(boolcheck->extended_table->elems[i], expected_table2[i]));
    }

    // Validate updated claim
    TEST_ASSERT_TRUE(f128_eq(boolcheck->claim, expected_eval));

    // Validate challenge is recorded
    TEST_ASSERT_EQUAL_UINT32(1, boolcheck->challenges->len);
    TEST_ASSERT_TRUE(f128_eq(boolcheck->challenges->elems[0], r));

    // Validate round_polys has one entry matching compressed_round_poly
    TEST_ASSERT_EQUAL_UINT32(1, boolcheck->round_polys_len);
    for (int i = 0; i < 3; ++i) {
        TEST_ASSERT_TRUE(f128_eq(boolcheck->round_polys->compressed_coeff[0].elems[i], compressed_round_poly->compressed_coeff->elems[i]));
    }

    // Manually set poly_coords to values from 0 to (4*128)
    size_t poly_coords_len = 4 * 128;
    if (boolcheck->poly_coords) free(boolcheck->poly_coords);
    boolcheck->poly_coords = points_init(poly_coords_len, f128_zero());
    // boolcheck->poly_coords->len = poly_coords_len;
    // boolcheck->poly_coords->elems = (F128*)malloc(sizeof(F128) * poly_coords_len);
    for (size_t i = 0; i < poly_coords_len; ++i) {
        boolcheck->poly_coords->elems[i] = f128_from_uint64(i);
    }

    // Call compute_algebraic(0, 1) and check output
    F128 alg_res[3];
    and_package_algebraic(&and_params, boolcheck->poly_coords, 0, 1, (F128*[]){alg_res, alg_res+1, alg_res+2});
    // These constants must match the Rust expected values
    F128 expected_alg[3] = {
        (F128){.low = 0x6ab015702ab00d03, .high = 0xe722c5702ab00070}, // 307232209640963015187300467897455873283
        (F128){.low = 0x35580ab815587601, .high = 0x929162b815580038}, // 194822172689813899653668252817418647041
        (F128){.low = 0x8002aaa2a002c2d7, .high = 0x32f512a2a002a022}  // 67733890487442795093175766325246739159
    };
    for (int i = 0; i < 3; ++i) {
        TEST_ASSERT_TRUE(f128_eq(alg_res[i], expected_alg[i]));
    }

    // Validate that the input polynomials are as expected
    for (int idx = 0; idx < 2; ++idx) {
        for (int j = 0; j < 8; ++j) {
            TEST_ASSERT_TRUE(f128_eq(polys->mle_poly[idx].coeffs[j], f128_from_uint64(j)));
        }
    }

    // Cleanup
    mle_poly_free(pq);
    mle_sequence_free(polys);
    points_free(points);
    bool_check_builder_free(builder);
    boolcheck_free(boolcheck);
    free(round_poly);
}