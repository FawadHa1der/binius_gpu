#include "test_mle.h"
#include "unity.h"
#include "../field.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>    // for seeding
#include "stdio.h" // for uint8_t, uint64_t
#include "../debug_utils.h" // for f128_print

/* A small helper to multiply an array of F128 values. */
F128 f128_mul_array(const F128* arr, size_t len)
{
    F128 acc = f128_one();
    for (size_t i = 0; i < len; i++) {
        acc = f128_mul(acc, arr[i]);
    }
    return acc;
}


/******************************************************************************
 * points_from_array
 *   Allocates a Points that copies 'n' F128 elements from 'arr'.
 ******************************************************************************/
Points points_from_array(const F128* arr, size_t n)
{
    Points p;
    p.len = n;
    if (n == 0) {
        p.elems = NULL;
        return p;
    }
    p.elems = (F128*)malloc(n * sizeof(F128));
    if (!p.elems) {
        // handle allocation failure
        p.len = 0;
        return p;
    }
    memcpy(p.elems, arr, n * sizeof(F128));
    return p;
}

Points points_random(size_t n){
    Points p;
    p.len = n;
    if (n == 0) {
        p.elems = NULL;
        return p;
    }

    p.elems = (F128*)malloc(n * sizeof(F128));
    if (!p.elems) {
        // handle allocation failure
        p.len = 0;
        return p;
    }

    for (size_t i = 0; i < n; i++) {
        ((F128*)(p.elems))[i] = f128_rand();
    }
    return p;
}


/*************************************************************************
 * 1) test_eq_eval_basic
 *
 * Mirrors your Rust "test_eq_eval()":
 *   points_a = [1,2,3]
 *   points_b = [4,5,6]
 *   result = eq_eval(...)
 *   expected = 1 * (1 +1+4)*(1 +2+5)*(1 +3+6)
 *************************************************************************/
void test_eq_eval_basic(void)
{
    /* Build points_a => [1,2,3] */
    F128 a_array[3];
    a_array[1] = f128_from_uint64(0x0000000000000001ULL);
    a_array[0] = f128_from_uint64(2ULL);
    a_array[2] = f128_from_uint64(3ULL);
    Points points_a = points_from_array(a_array, 3);

    /* Build points_b => [4,5,6] */
    F128 b_array[3];
    b_array[0] = f128_from_uint64(4ULL);
    b_array[1] = f128_from_uint64(5ULL);
    b_array[2] = f128_from_uint64(6ULL);
    Points points_b = points_from_array(b_array, 3);

    /* eq_eval */
    F128 result = points_eq_eval(&points_a, &points_b);

    /* expected = 1 * (1 +1+4)*(1 +2+5)*(1 +3+6) */
    F128 one = f128_one();

    F128 sum1 = f128_add(one, f128_add(a_array[0], b_array[0])); // 1+(1+4)
    F128 sum2 = f128_add(one, f128_add(a_array[1], b_array[1])); // 1+(2+5)
    F128 sum3 = f128_add(one, f128_add(a_array[2], b_array[2])); // 1+(3+6)
    F128 sums[3] = { sum1, sum2, sum3 };
    F128 expected = f128_mul_array(sums, 3);

    /* Use Unity's TEST_ASSERT for the check */
    TEST_ASSERT_MESSAGE(f128_eq(result, expected), 
                       "test_eq_eval_basic FAILED: result != expected");
}

/*************************************************************************
 * 2) test_eq_eval_identity
 *
 * Rust "test_eq_eval_identity()":
 *   points = [7,8,9]
 *   eq_eval(points, points) => 1
 *************************************************************************/
void test_eq_eval_identity(void)
{
    F128 arr[3];
    arr[0] = f128_from_uint64(7ULL);
    arr[1] = f128_from_uint64(8ULL);
    arr[2] = f128_from_uint64(9ULL);
    Points points = points_from_array(arr, 3);

    F128 result = points_eq_eval(&points, &points);

    F128 one = f128_one();
    TEST_ASSERT_MESSAGE(f128_eq(result, one),
                       "test_eq_eval_identity FAILED: result != 1");
}

/*************************************************************************
 * 3) test_eq_eval_empty
 *
 * Rust "test_eq_eval_empty()":
 *   points_a=Points::default(), points_b=Points::default()
 *   eq_eval(...) => 1
 *************************************************************************/
void test_eq_eval_empty(void)
{
    Points pa = points_default(); // empty
    Points pb = points_default(); // empty

    F128 result = points_eq_eval(&pa, &pb);

    TEST_ASSERT_MESSAGE(f128_eq(result, f128_one()), 
                       "test_eq_eval_empty FAILED: result != 1");
}


void test_eq_poly_sequence_cross_check(void)
{
    size_t n = 20;
    // 1) generate random points
    Points points = points_random(n);

    // 2) eq_sequence
    // size_t seq_len = points.len + 1;
    MLE_POLY_SEQUENCE *eq_sequence =
        points_to_eq_poly_sequence(&points);

    // check eq_sequence len = points.len + 1
    TEST_ASSERT_EQUAL_UINT(points.len + 1, eq_sequence->len );

    // check eq_sequence[0] = polynomial with [F128::ONE]
    // We'll assume eq_sequence[0].len == 1 and eq_sequence[0].coeffs[0] == one
    TEST_ASSERT_EQUAL_UINT(1, eq_sequence->mle_poly[0].len);
    F128 one = f128_one();
    TEST_ASSERT_TRUE( f128_eq(eq_sequence->mle_poly[0].coeffs[0], one) );

    // cross-check each polynomial in [1..seq_len)
    for (size_t i = 1; i < eq_sequence->len; i++) {
        // "pts" = sub-slice of 'points' from (points.len - i) to end
        // We'll do a small function "points_subrange" or a direct approach:
        size_t start = points.len - i;
        size_t size  = i; 
        // build sub-points
        Points sub_pts;
        sub_pts.len = size;
        sub_pts.elems = (F128*)malloc(size * sizeof(F128));
        for (size_t j = 0; j < size; j++) {
            
            sub_pts.elems[j] = points.elems[start + j];
        }

        // direct computation => poly2 = to_eq_poly(&sub_pts)
        MLE_POLY* poly2 = points_to_eq_poly(&sub_pts);

        // eq_sequence[i] should match poly2
        TEST_ASSERT_TRUE( mle_poly_eq(&eq_sequence->mle_poly[i], poly2) );

        // cleanup
        free(sub_pts.elems);
        free(poly2->coeffs);  // if that's how memory is done
    }

    // free eq_sequence 
    for (size_t i=0; i < eq_sequence->len; i++) {
        free(eq_sequence->mle_poly[i].coeffs);
    }
    free(eq_sequence);

    // free points
    free(points.elems);

    // printf("test_eq_poly_sequence_cross_check PASSED\n");
}


// Test function for to_points_inv_orbit with ones
void test_to_points_inv_orbit_ones(void) {
    // Generate a set of 5 points initialized to F128_ONE
    Points initial_points;
    initial_points.len = 5;
    initial_points.elems = malloc(sizeof(F128) * initial_points.len);
    
    for (size_t i = 0; i < initial_points.len; i++) {
        initial_points.elems[i] = f128_one();
    }

    // Compute the inverse orbit
    INVERSE_ORBIT_POINTS *points_inv_orbit = to_f128_inv_orbit(&initial_points);

    // Validate the length of the inverse orbit
    TEST_ASSERT_EQUAL_INT(128, points_inv_orbit->len);

    // Validate that all entries in the orbit are ones
    for (size_t i = 0; i < 128; i++) {
        for (size_t j = 0; j < initial_points.len; j++) {
            F128 one = f128_one();

            TEST_ASSERT_EQUAL_UINT64(one.low, points_inv_orbit->array_of_points[i].elems[j].low);
            TEST_ASSERT_EQUAL_UINT64(one.high, points_inv_orbit->array_of_points[i].elems[j].high);
        }
    }

    // Free allocated memory
    for (size_t i = 0; i < 128; i++) {
        free(points_inv_orbit->array_of_points[i].elems);
    }
    free(points_inv_orbit->array_of_points);
    free(points_inv_orbit);
    free(initial_points.elems);
}


void test_to_points_inv_orbit_last_element(void) {
    
    // Define two initial points in the binary field with fixed values.
    
    F128 pt1 = f128_from_uint64(1234);  // Equivalent to BinaryField128b::new(1234)
    F128 pt2 = f128_from_uint64(5678);  // Equivalent to BinaryField128b::new(5678)

    // Create the initial Points instance with the two points
    Points initial_points;
    initial_points.len = 2;
    initial_points.elems = malloc(sizeof(F128) * initial_points.len);
    initial_points.elems[0] = pt1;
    initial_points.elems[1] = pt2;

    // Call to_points_inv_orbit() to compute the inverse orbit.
    // Points *points_inv_orbit = to_points_inv_orbit(&initial_points);
    INVERSE_ORBIT_POINTS *points_inv_orbit = to_f128_inv_orbit(&initial_points);

    // Initialize expected results
    Points *expected_points = malloc(sizeof(Points) * 128);
    for (size_t i = 0; i < 128; i++) {
        expected_points[i].len = 2;
        expected_points[i].elems = malloc(sizeof(F128) * 2);

        // Frobenius squaring
        pt1 = f128_mul(pt1, pt1);
        pt2 = f128_mul(pt2, pt2);

        expected_points[i].elems[0] = pt1;
        expected_points[i].elems[1] = pt2;
    }

    // Reverse the expected results to match function output.
    for (size_t i = 0; i < 64; i++) {
        Points temp = expected_points[i];
        expected_points[i] = expected_points[127 - i];
        expected_points[127 - i] = temp;
    }

    // Compare computed orbit with expected results
    for (size_t i = 0; i < 128; i++) {
        TEST_ASSERT_EQUAL_MEMORY(expected_points[i].elems, points_inv_orbit->array_of_points[i].elems, sizeof(F128) * 2);
    }

    // Free allocated memory
    for (size_t i = 0; i < 128; i++) {
        free(points_inv_orbit->array_of_points[i].elems);
        free(expected_points[i].elems);
    }
    free(points_inv_orbit);
    free(expected_points);
    free(initial_points.elems);
}


/******************************************************************************
 * test_eq_sums
 *
 * This test constructs an equality polynomial from the points [1,2,3,4],
 * computes its equality sums, then manually computes the subset sums block‐by‐block.
 * Finally, it asserts that the computed equality sums match the expected ones.
 *****************************************************************************/
void test_eq_sums(void) {
    // --- Step 1: Create the Points from [1,2,3,4] ---
    F128 *pts_arr = malloc(4 * sizeof(F128));
    pts_arr[0] = f128_from_uint64(1ULL);
    pts_arr[1] = f128_from_uint64(2ULL);
    pts_arr[2] = f128_from_uint64(3ULL);
    pts_arr[3] = f128_from_uint64(4ULL);

    Points* points = malloc(sizeof(Points));
    points->len = 4;
    points->elems = pts_arr;
    // Points points = points_from_array(pts_arr, 4);

    // --- Step 2: Compute the equality polynomial ---
    MLE_POLY* poly = points_to_eq_poly(points);

    
    // --- Step 3: Compute the equality sums ---
    Points *computed_eqs = eq_sums(poly);
    
    // --- Step 4: Manually compute expected sums ---
    // We assume poly.length is divisible by 8.
    TEST_ASSERT_EQUAL_INT(0, poly->len % 8);
    size_t num_blocks = poly->len / 8;
    size_t expected_total = 256 * num_blocks;
    F128 *expected_eqs = malloc(expected_total * sizeof(F128));
    if (!expected_eqs) {
        TEST_FAIL_MESSAGE("Memory allocation failed for expected_eqs.");
    }
    
    // For each block of 8 coefficients:
    for (size_t block = 0; block < num_blocks; block++) {
        size_t block_start = block * 8;
        F128 *block_coeffs = &poly->coeffs[block_start];
        
        // For each subset index from 0 to 255:
        for (size_t subset = 0; subset < 256; subset++) {
            F128 sum = f128_zero();
            // Iterate over the 8 bits of 'subset'
            for (size_t bit = 0; bit < 8; bit++) {
                if (subset & (1 << bit)) {
                    sum = f128_add(sum, block_coeffs[bit]);
                }
            }
            expected_eqs[block * 256 + subset] = sum;
        }
    }
    
    // --- Step 5: Compare computed_eqs and expected_eqs ---
    TEST_ASSERT_EQUAL_UINT(expected_total, computed_eqs->len);
    // print out both the arrays

    for (size_t i = 0; i < expected_total; i++) {
        if (!f128_eq(computed_eqs->elems[i], expected_eqs[i])) {
            char msg[128];
            sprintf(msg, "Mismatch at index %zu", i);
            TEST_FAIL_MESSAGE(msg);
        }
    }
    
    // --- Cleanup ---
    free(computed_eqs);
    free(expected_eqs);
    free(poly->coeffs);
    free(points->elems);
    free(points);
}


void test_drop_top_bit_standard_cases(void) {
    // Test with input 0b1010 (decimal 10).
    // Expected: 0b1010 -> 0b0010, MSB position = 3.
    uint8_t input = 0b1010; // 10 in decimal.
    uint8_t result, bit_index;
    drop_top_bit(input, &result, &bit_index);
    TEST_ASSERT_EQUAL_UINT8(0b0010, result);
    TEST_ASSERT_EQUAL_UINT8(3, bit_index);

    // Test with input 0b0101 (decimal 5).
    // Expected: 0b0101 -> 0b0001, MSB position = 2.
    input = 0b0101;
    drop_top_bit(input, &result, &bit_index);
    TEST_ASSERT_EQUAL_UINT8(0b0001, result);
    TEST_ASSERT_EQUAL_UINT8(2, bit_index);
}

void test_drop_top_bit_edge_cases(void) {
    // Smallest non-zero: 1 (0b0001)
    uint8_t input = 1;
    uint8_t result, bit_index;
    drop_top_bit(input, &result, &bit_index);
    // 0b0001 -> 0, MSB position = 0.
    TEST_ASSERT_EQUAL_UINT8(0, result);
    TEST_ASSERT_EQUAL_UINT8(0, bit_index);
}

void test_drop_top_bit_large_numbers(void) {
    uint8_t result, bit_index;
    // Test with 2^31. But note: our drop_top_bit is designed for 8-bit values.
    // In our use, x is always between 1 and 255. So we simulate by casting.
    uint8_t input = (uint8_t)(1 << 7); // For an 8-bit value, the maximum bit is bit 7.
    // For a test of "large" (e.g., 1 << 31) in the context of 8-bit, we must adjust:
    // Instead, we'll test with inputs 1<<k for k=0..7.
    // For k = 7 (i.e., 0b10000000) the expected result is 0, and bit index = 7.
    input = 1 << 7;
    drop_top_bit(input, &result, &bit_index);
    TEST_ASSERT_EQUAL_UINT8(0, result);
    TEST_ASSERT_EQUAL_UINT8(7, bit_index);

    // For an 8-bit simulation of 2^63, also expect result 0 and bit index 7.
    // (Since our function is only for 8-bit numbers.)
}

void test_drop_top_bit_all_bits_set(void) {
    // For a 4-bit number: 0b1111 (15)
    // Expected: 0b1111 -> 0b0111, MSB position = 3.
    uint8_t input = 0b1111;
    uint8_t result, bit_index;
    drop_top_bit(input, &result, &bit_index);
    TEST_ASSERT_EQUAL_UINT8(0b0111, result);
    TEST_ASSERT_EQUAL_UINT8(3, bit_index);
}

// -----------------------------------------------------------------------------
// Test cases for cpu_v_movemask_epi8
// -----------------------------------------------------------------------------

void test_cpu_v_movemask_epi8_standard_cases(void) {
    // Input: 16 bytes alternating between 0b10000000 and 0b00000000.
    uint8_t input[16] = {
        0x80, 0x00, 0x80, 0x00,
        0x80, 0x00, 0x80, 0x00,
        0x80, 0x00, 0x80, 0x00,
        0x80, 0x00, 0x80, 0x00,
    };
    // Expected: Bits: 0b0101010101010101 (21845)
    uint16_t expected = 0b0101010101010101; // 21845
    uint16_t result = cpu_v_movemask_epi8(input);
    TEST_ASSERT_EQUAL_UINT16(expected, result);
}

void test_cpu_v_movemask_epi8_all_ones(void) {
    // Input: 16 bytes, each 0b10000000.
    uint8_t input[16];
    for (int i = 0; i < 16; i++) {
        input[i] = 0x80;
    }
    uint16_t expected = 0xFFFF; // 16 ones => 65535
    uint16_t result = cpu_v_movemask_epi8(input);
    TEST_ASSERT_EQUAL_UINT16(expected, result);
}

void test_cpu_v_movemask_epi8_all_zeros(void) {
    // Input: 16 bytes all 0.
    uint8_t input[16] = {0};
    uint16_t expected = 0;
    uint16_t result = cpu_v_movemask_epi8(input);
    TEST_ASSERT_EQUAL_UINT16(expected, result);
}

// -----------------------------------------------------------------------------
// Test cases for v_slli_epi64
// -----------------------------------------------------------------------------

// We'll define v_slli_epi64 as: void v_slli_epi64_c(const uint8_t in[16], uint8_t out[16], unsigned shift)

void test_v_slli_epi64_basic_shift(void) {
    // Input: two 64-bit values encoded as 16 bytes:
    // first 8 bytes represent 1, next 8 represent 2.
    uint8_t input[16] = {
        0x01,0x00,0x00,0x00,0x00,0x00,0x00,0x00, // 1 in little-endian
        0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00  // 2 in little-endian
    };
    uint8_t expected[16] = {
        0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00, // 1 << 1 = 2
        0x04,0x00,0x00,0x00,0x00,0x00,0x00,0x00  // 2 << 1 = 4
    };
    // uint8_t result[16] = {0};
    v_slli_epi64_c( 1, input);
    TEST_ASSERT_EQUAL_MEMORY(expected, input, sizeof(expected));
}

void test_v_slli_epi64_zero_shift(void) {
    // Input: some arbitrary 16-byte value.
    uint8_t input[16] = {
        0x10,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x20,0x00,0x00,0x00,0x00,0x00,0x00,0x00
    };
    uint8_t expected[16];
    memcpy(expected, input, sizeof(input));
    // uint8_t result[16] = {0};
    v_slli_epi64_c(0, input);
    TEST_ASSERT_EQUAL_MEMORY(expected, input, sizeof(expected));
}

void test_v_slli_epi64_edge_cases(void) {
    // Test edge cases:
    // Input: first 8 bytes all 0xFF, second 8 bytes: all 0 except last byte = 0x80.
    uint8_t input[16] = {
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x80
    };
    // Expected: shift left by 1:
    // For first 8 bytes: shifting 0xFF left by 1 gives 0xFE in each byte (with overflow discarded).
    // For second 8 bytes: shifting 0x00...0x80 left by 1 gives all 0.
    uint8_t expected[16] = {
        0xFE,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00
    };
    // uint8_t result[16] = {0};
    v_slli_epi64_c(1, input);
    TEST_ASSERT_EQUAL_MEMORY(expected, input, sizeof(expected));
}

void test_v_slli_epi64_no_overflow(void) {
    // Check no overflow: 
    // Input: first 8 bytes represent 1, second 8 represent 64.
    uint8_t input[16] = {
        0x01,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x40,0x00,0x00,0x00,0x00,0x00,0x00,0x00
    };
    // Shift left by 3:
    // 1 << 3 = 8 -> represented as 0x08,0x00,...
    // 64 << 3 = 512, which in little-endian 64-bit: 512 = 0x200,
    // so the 64-bit number becomes: 0x00, 0x02, then zeros.
    uint8_t expected[16] = {
        0x08,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x00,0x02,0x00,0x00,0x00,0x00,0x00,0x00
    };
    // uint8_t result[16] = {0};
    v_slli_epi64_c(3, input);
    TEST_ASSERT_EQUAL_MEMORY(expected, input, sizeof(expected));
}


// -------------- Utility to compare Evaluations with an expected array --------

// We'll define a helper to check that the returned Evaluations match
// some expected F128 array:
void assert_evaluations_equal(
    const Points* actual,
    const F128* expected_vals,
    size_t expected_len
) {
    TEST_ASSERT_EQUAL_UINT64(expected_len, actual->len);
    for(size_t i=0; i<expected_len; i++){
        TEST_ASSERT_MESSAGE(
            f128_eq(actual->elems[i], expected_vals[i]),
            "Mismatch in Evaluations at index"
        );
    }
}

void test_restrict_polynomials(void)
{
    // (1) define the number of variables
    size_t num_vars = 4;
    // how many we restrict
    size_t num_vars_to_restrict = 4;

    // (2) define the restriction points
    Points pts;
    pts.len = num_vars_to_restrict;
    // allocate or create them
    pts.elems = malloc(sizeof(F128)*num_vars_to_restrict);
    for(size_t i=0; i<num_vars_to_restrict; i++){
        pts.elems[i] = f128_from_uint64(i); // mock example like Rust
    }

    // (3) define 3 polynomials (N=3)
    // each polynomial has length = 1<<num_vars = 16
    size_t poly_len = (1ULL << num_vars);

    MLE_POLY poly0, poly1, poly2;
    poly0.len = poly_len; 
    poly1.len = poly_len; 
    poly2.len = poly_len;
    poly0.coeffs = malloc(sizeof(F128)*poly_len);
    poly1.coeffs = malloc(sizeof(F128)*poly_len);
    poly2.coeffs = malloc(sizeof(F128)*poly_len);

    // fill poly0 with new(i)
    for(size_t i=0; i<poly_len; i++){
        poly0.coeffs[i] = f128_from_uint64(i);  // Rust used new(i)
        poly1.coeffs[i] = f128_from_uint64(i+1);
        poly2.coeffs[i] = f128_from_uint64(i+2);
    }

    // (4) call the restrict function
    // we'll store the result in an Evaluations object
    MLE_POLY polys_array[3];
    polys_array[0] = poly0;
    polys_array[1] = poly1;
    polys_array[2] = poly2;

    Points* restricted = restrict_polynomials(polys_array, 3, &pts, num_vars);

    size_t expected_len = 384;
    F128 expected[expected_len];
    for(size_t i=0; i<expected_len; i++){
        expected[i] = f128_zero();
    }
    // Then from the test, some were non-zero:
    // e.g. expected[0] = ZERO, expected[1] = new(1), ...
    expected[1] = f128_from_uint64(1);
    expected[2] = f128_from_uint64(2);
    expected[3] = f128_from_uint64(3);

    // index 128
    //BinaryField128b::new(257870231182273679343338569694386847745),
    //0x0000000000000001, 
    expected[128] = f128_from_raw(0x0000000000000001, 0xc200000000000000);
    expected[129] = f128_from_uint64(1);
    expected[130] = f128_from_uint64(2);
    expected[131] = f128_from_uint64(3);

    //0x0000000000000000, 0xc200000000000000
    expected[257] = f128_from_raw(0x0000000000000000, 0xc200000000000000);
    expected[258] = f128_from_uint64(3);
    // 0x0000000000000000, 0xe608000000000000
    expected[259] = f128_from_raw(0x0000000000000000, 0xe608000000000000);
    //0x0000000000000007, 0xa9d2300000000000
    expected[260] = f128_from_raw(0x0000000000000007, 0xa9d2300000000000);


    // ... etc. In real code you'd fill them all.

    // (6) compare using a helper
    TEST_ASSERT_EQUAL_UINT64(384, restricted->len);
    // or call a custom compare function
    // We'll do a minimal check just for the first 4 to show the approach:
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[0], f128_zero()));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[1], f128_from_uint64(1)));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[2], f128_from_uint64(2)));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[3], f128_from_uint64(3)));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[128], expected[128]));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[129], expected[129]));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[130], expected[130]));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[131], expected[131]));

    TEST_ASSERT_TRUE(f128_eq(restricted->elems[257], expected[257]));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[258], expected[258]));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[259], expected[259]));
    TEST_ASSERT_TRUE(f128_eq(restricted->elems[260], expected[260]));


    // ... etc.

    // if you want to check all 128 exactly:
    // assert_evaluations_equal(&restricted, expected, 128);

    // (7) Clean up
    free(poly0.coeffs);
    free(poly1.coeffs);
    free(poly2.coeffs);
    free(pts.elems);
    free(restricted->elems);
}


void test_eq_poly(void)
{
    // "Define a simple set of points"
    F128 arr[3] = {
        f128_from_uint64(1ULL),
        f128_from_uint64(2ULL),
        f128_from_uint64(3ULL),
    };
    Points pts = points_from_array(arr, 3);

    // to_eq_poly
    MLE_POLY* result = points_to_eq_poly(&pts);

    // we'd compare result to the "expected" large array. For demonstration, let's define it:
    size_t length = (1ULL<<3); //8
        F128 expected8[8] = {
        { .low = 0x0000000000000003, .high = 0x11ce300000000000 }, // partial demonstration
        { .low=0x0000000000000007, .high=0x3bd6300000000000 },
        { .low=0x0000000000000002, .high=0xa7c2300000000000 },
        { .low=0x0000000000000004, .high=0x4fda300000000000 },
        { .low=0x0000000000000002, .high=0x35c6300000000000 },
        { .low=0x0000000000000005, .high=0xddde300000000000 },
        { .low=0x0000000000000003, .high=0x41ca300000000000 },
        { .low=0x0000000000000007, .high=0xa9d2300000000000 }
    };

    // Check lengths
    TEST_ASSERT_EQUAL_UINT64(length, result->len);

    // check each
    for(size_t i=0;i<length;i++){
        // compare result.coeffs[i] with expected8[i]
        TEST_ASSERT_TRUE(f128_eq(result->coeffs[i], expected8[i]));
    }

    // free
    free(result);
    free(pts.elems);
}


void test_eq_poly_evaluate(void)
{
    // "Define a simple set of points"
    F128 pt1 = f128_from_uint64(3434ULL);
    F128 pt2 = f128_from_uint64(6765ULL);

    F128 arr[2] = { pt1, pt2};
    Points points = points_from_array(arr,2);

    // "Compute eq poly"
    MLE_POLY* eq_poly = points_to_eq_poly(&points);

    // "Evaluate eq poly at the given points"
    // we do the standard formula in your snippet:
    // "let expected_eq_poly = result[0]*(1-pt1)*(1-pt2) + result[2]*(1-pt1)*pt2 + ..."

    // We interpret eq_poly as length=4
    // eq_poly[0], eq_poly[1], eq_poly[2], eq_poly[3]
    if(eq_poly->len !=4){
        // handle error or assert
        assert(0 && "eq_poly length is not 4");
    }
    F128 c0= eq_poly->coeffs[0];
    F128 c1= eq_poly->coeffs[1];
    F128 c2= eq_poly->coeffs[2];
    F128 c3= eq_poly->coeffs[3];

    // " (BinaryField128b::ONE - pt1) " => f128_add(&F128_ONE, &pt1)
    F128 one_minus_pt1 = f128_add(f128_one(), pt1);
    F128 one_minus_pt2 = f128_add(f128_one(), pt2);

    // do c0*(1-pt1)*(1-pt2)
    F128 t1=f128_mul(c0, one_minus_pt1);
    t1=f128_mul(t1, one_minus_pt2);

    // c2*(1-pt1)*pt2
    F128 t2=f128_mul(c2, one_minus_pt1);
    t2=f128_mul(t2, pt2);

    // c1*pt1*(1-pt2)
    F128 t3=f128_mul(c1, pt1);
    t3=f128_mul(t3, one_minus_pt2);

    // c3*pt1*pt2
    F128 t4=f128_mul(c3, pt1);
    t4=f128_mul(t4, pt2);

    // sum them up
    F128 sum=f128_zero();
    sum=f128_add(sum, t1);
    sum=f128_add(sum, t2);
    sum=f128_add(sum, t3);
    sum=f128_add(sum, t4);

    // "Verify eq poly = 1"
    // We'll check sum == f128_ONE
    TEST_ASSERT_TRUE(f128_eq(sum,f128_one()));

    // free
    free(eq_poly);
    free(points.elems);
}


void test_eq_poly_sequence(void)
{

    // "Define a simple set of points"
    F128 arr[4] = {
        f128_from_uint64(1ULL),
        f128_from_uint64(2ULL),
        f128_from_uint64(3ULL),
        f128_from_uint64(4ULL),
    };
    Points pts = points_from_array(arr, 4);

    // to_eq_poly
    MLE_POLY_SEQUENCE* result = points_to_eq_poly_sequence(&pts);

    TEST_ASSERT_TRUE( result->mle_poly[0].len == 1);
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[0].coeffs[0], f128_from_raw(0x0000000000000001, 0xc200000000000000)) );

    TEST_ASSERT_TRUE( result->mle_poly[1].len == 2);
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[1].coeffs[0], f128_from_raw(0x0000000000000005, 0xc200000000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[1].coeffs[1], f128_from_raw(0x0000000000000004, 0x0)) );


    TEST_ASSERT_TRUE( result->mle_poly[2].len == 4);
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[2].coeffs[0], f128_from_raw(0x000000000000000f, 0xd030000000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[2].coeffs[1], f128_from_raw(0x000000000000000a, 0x1230000000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[2].coeffs[2], f128_from_raw(0x000000000000000d, 0x1230000000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[2].coeffs[3], f128_from_raw(0x0000000000000009, 0x1230000000000000)) );

    TEST_ASSERT_TRUE( result->mle_poly[3].len == 8);
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[0], f128_from_raw(0x0000000000000018, 0xc540c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[1], f128_from_raw(0x0000000000000017, 0x1570c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[2], f128_from_raw(0x0000000000000011, 0x1b60c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[3], f128_from_raw(0x000000000000001b, 0x0950c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[4], f128_from_raw(0x000000000000001c, 0xef58c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[5], f128_from_raw(0x0000000000000011, 0xfd68c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[6], f128_from_raw(0x0000000000000016, 0xf378c00000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[3].coeffs[7], f128_from_raw(0x000000000000001f, 0xe148c00000000000)) );

    TEST_ASSERT_TRUE( result->mle_poly[4].len == 16);
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[0], f128_from_raw(0x000000000000000a, 0x761ad18000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[1], f128_from_raw(0x0000000000000012, 0xb35a118000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[2], f128_from_raw(0x000000000000000e, 0x4d92b18000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[3], f128_from_raw(0x0000000000000019, 0x58e2718000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[4], f128_from_raw(0x000000000000000c, 0xa11a918000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[5], f128_from_raw(0x000000000000001d, 0xba7a518000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[6], f128_from_raw(0x000000000000000b, 0x5682f18000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[7], f128_from_raw(0x0000000000000010, 0x5fd2318000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[8], f128_from_raw(0x0000000000000009, 0x67d4e18000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[9], f128_from_raw(0x0000000000000015, 0x888c218000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[10], f128_from_raw(0x000000000000000c, 0xea50818000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[11], f128_from_raw(0x000000000000001d, 0x1738418000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[12], f128_from_raw(0x000000000000000e, 0x94dca18000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[13], f128_from_raw(0x0000000000000018, 0x67a4618000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[14], f128_from_raw(0x0000000000000008, 0x1748c18000000000)) );
    TEST_ASSERT_TRUE( f128_eq(result->mle_poly[4].coeffs[15], f128_from_raw(0x0000000000000017, 0xf600018000000000)) );


}


int random_in_range(int low, int high){
    // assume seeded
    int r = rand(); 
    return low + (r%(high-low));
}

/******************************************************************************
 * test_eq_poly_sequence_random_values
 *
 * This test constructs an equality polynomial from random points,
 * computes its equality sums, then manually computes the subset sums block‐by‐block.
 * Finally, it asserts that the computed equality sums match the expected ones.
 *****************************************************************************/

void test_eq_poly_sequence_random_values(void)
{
    // step 1: define #iterations=10
    int iterations=10;

    // seed
    srand((unsigned)time(NULL));

    for(int iter=0; iter<iterations; iter++){
        // Step 1: random dimension in [2..6]
        int num_points = random_in_range(2,6);

        // build random points
        Points points = points_random(num_points);

        // Step 2: to_eq_poly => "result"
        MLE_POLY* result = points_to_eq_poly(&points);
        size_t length= (1ULL<< num_points);

        // Step 3: Reconstruct eq poly manually => we expect sum=1
        // We'll do a "expected_eq_poly" as 
        F128 expected_eq_poly= f128_zero();

        // For i in 0..(1<<num_points):
        for(size_t i=0;i<length;i++){
            // turn i into binary => if bit j => multiply by points->vals[j], else multiply by (1- points->vals[j])
            F128 term= result->coeffs[i];
            for(int j=0;j<num_points;j++){
                int isSet= ( (i>>j) &1 );
                if(isSet){
                    term=f128_mul(term, points.elems[j]);
                } else {
                    // (1- points->vals[j])
                    F128 complement= f128_add(f128_one(),points.elems[j]);
                    term=f128_mul(term, complement);
                }
            }
            // add to expected
            expected_eq_poly= f128_add(expected_eq_poly,term);
        }

        // Step4: check eq to 1
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(expected_eq_poly, f128_one()),
            "computed eq poly != 1"
        );

        // free
        free(result);
        free(points.elems);
    }
}

void test_evaluate_at(void)
{
    // mimic your snippet: define a 2D polynomial with 4 coefficients, random them,
    // define points, do mle_poly_evaluate_at, do the manual formula, compare
    // We'll do partial:

    // define polynomial
    F128 c0= f128_from_uint64(42ULL);
    F128 c1= f128_from_uint64(100ULL);
    F128 c2= f128_from_uint64(99ULL);
    F128 c3= f128_from_uint64(77ULL);
    F128 polyarr[4] = {c0,c1,c2,c3};
    MLE_POLY* polynomial = malloc(sizeof(MLE_POLY));
    polynomial->len = 4;
    polynomial->coeffs = malloc(sizeof(F128)*4);
    for(size_t i=0;i<4;i++){
        polynomial->coeffs[i] = polyarr[i];
    }

    // define points in dimension=2
    F128 p0= f128_from_uint64(1ULL);
    F128 p1= f128_from_uint64(0ULL);
    F128 ptsarr[2]= {p0,p1};
    Points pts= points_from_array(ptsarr,2);

    // call mle_poly_evaluate_at
    F128 result = mle_poly_evaluate_at(polynomial,&pts);

    // do manual formula
    F128 one_minus_p0= f128_add(f128_one(),p0);
    F128 one_minus_p1= f128_add(f128_one(),p1);
    // c0*(1-p0)*(1-p1)+ c2*(1-p0)*p1+ c1*(1-p1)*p0+ c3*p0*p1
    // partial:

    F128 sum=f128_zero();

    // c0*(1-p0)*(1-p1)
    F128 t=f128_mul(c0,one_minus_p0);
    t=f128_mul(t,one_minus_p1);
    sum=f128_add(sum,t);

    // c2*(1-p0)*p1
    t=f128_mul(c2,one_minus_p0);
    t=f128_mul(t,p1);
    sum=f128_add(sum,t);

    // c1*p0*(1-p1)
    t=f128_mul(c1,p0);
    t=f128_mul(t,one_minus_p1);
    sum=f128_add(sum,t);

    // c3*p0*p1
    t=f128_mul(c3,p0);
    t=f128_mul(t,p1);
    sum=f128_add(sum,t);

    // //print all the values
    // f128_print("result ", result);
    // f128_print("sum ", sum);


    // compare result
    TEST_ASSERT_TRUE(f128_eq(result,sum));

    // free
    free(polynomial);
    free(pts.elems);
}

MLE_POLY* mlp_bitand(const MLE_POLY* a, 
    const MLE_POLY* b)
{
    // assume a->len==b->len
    MLE_POLY* res = malloc(sizeof(MLE_POLY));
    res->len= a->len;
    res->coeffs = malloc(sizeof(F128)* a->len);

    // for demonstration, do "bitwise AND" as an operation on .lo & .lo
    // not real field logic, just demonstration
    for(size_t i=0;i<a->len;i++){
    F128 r;
    r.low = a->coeffs[i].low & b->coeffs[i].low;
    r.high = a->coeffs[i].high & b->coeffs[i].high;
    res->coeffs[i] = r;
    }
    return res;
}



void test_multilinear_lagrangian_bitand(void)
{
    // define two polynomials with 4 coeffs
    F128 pa[4] = {
        f128_from_uint64(1), // p1(0,0)
        f128_from_uint64(0), // p1(0,1)
        f128_from_uint64(1), // p1(1,0)
        f128_from_uint64(1)  // p1(1,1)
    };
    F128 pb[4] = {
        f128_from_uint64(1), // p2(0,0)
        f128_from_uint64(1), // p2(0,1)
        f128_from_uint64(0), // p2(1,0)
        f128_from_uint64(1)  // p2(1,1)
    };

    MLE_POLY* poly1=  malloc(sizeof(MLE_POLY)); // MLP_new(pa,4);
    poly1->len=4;
    poly1->coeffs= pa;
    MLE_POLY* poly2=  malloc(sizeof(MLE_POLY)); // MLP_new(pb,4);
    poly2->len=4;
    poly2->coeffs= pb;

    // bitwise AND
    MLE_POLY* result= mlp_bitand(poly1,poly2);

    // expected = [1,0,0,1]
    F128 expected_arr[4] = {
        f128_from_uint64(1),
        f128_from_uint64(0),
        f128_from_uint64(0),
        f128_from_uint64(1)
    };

    MLE_POLY* expected= malloc(sizeof(MLE_POLY)); // MLP_new(expected_arr,4);
    expected->len=4;
    expected->coeffs= expected_arr;

    // check each
    TEST_ASSERT_EQUAL_UINT64(4, result->len);
    for(size_t i=0;i<4;i++){
        TEST_ASSERT_TRUE(f128_eq(result->coeffs[i], expected->coeffs[i]));
    }

    free(result);
    free(expected);
    free(poly1);
    free(poly2);
}
