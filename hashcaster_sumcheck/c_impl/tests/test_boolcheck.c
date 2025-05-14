#include "test_boolcheck.h"
#include "unity.h"


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
    F128* result = extend_n_tables(polys, N, dims, C, trit_map, builder,
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
        TEST_ASSERT_TRUE(f128_eq(result[i], expected[i]));
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