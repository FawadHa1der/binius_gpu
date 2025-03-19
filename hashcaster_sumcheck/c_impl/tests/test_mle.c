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
        p.elems[i] = f128_rand();
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

    printf("test_eq_poly_sequence_cross_check PASSED\n");
}
