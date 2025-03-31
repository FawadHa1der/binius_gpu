#include "test_mle.h"
#include "unity.h"
#include "../field.h"
#include "stdlib.h"
#include "string.h"


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
        to_eq_poly_sequence(&points);

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
