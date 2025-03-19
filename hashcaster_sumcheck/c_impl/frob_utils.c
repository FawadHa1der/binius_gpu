#include "frob_utils.h"
#include "field.h"

static void test_eq_poly_sequence_cross_check(void)
{
    // We'll assume we have a random generator, or stub. Pass NULL if we don't:
    void *rng = NULL; 
    size_t n = 20;
    // 1) generate random points
    Points points = points_random(n, rng);

    // 2) eq_sequence
    size_t seq_len = 0;
    MultilinearLagrangianPolynomial *eq_sequence =
        to_eq_poly_sequence(&points, &seq_len);

    // check eq_sequence length = points.length + 1
    TEST_ASSERT_EQUAL_UINT(points.length + 1, seq_len);

    // check eq_sequence[0] = polynomial with [F128::ONE]
    // We'll assume eq_sequence[0].length == 1 and eq_sequence[0].coeffs[0] == one
    TEST_ASSERT_EQUAL_UINT(1, eq_sequence[0].length);
    F128 one = f128_one();
    TEST_ASSERT_TRUE( f128_eq(eq_sequence[0].coeffs[0], one) );

    // cross-check each polynomial in [1..seq_len)
    for (size_t i = 1; i < seq_len; i++) {
        // "pts" = sub-slice of 'points' from (points.length - i) to end
        // We'll do a small function "points_subrange" or a direct approach:
        size_t start = points.length - i;
        size_t size  = i; 
        // build sub-points
        Points sub_pts;
        sub_pts.length = size;
        sub_pts.data = (F128*)malloc(size * sizeof(F128));
        for (size_t j = 0; j < size; j++) {
            sub_pts.data[j] = points.data[start + j];
        }

        // direct computation => poly2 = to_eq_poly(&sub_pts)
        MultilinearLagrangianPolynomial poly2 = to_eq_poly(&sub_pts);

        // eq_sequence[i] should match poly2
        TEST_ASSERT_TRUE( mlp_eq(&eq_sequence[i], &poly2) );

        // cleanup
        free(sub_pts.data);
        free(poly2.coeffs);  // if that's how memory is done
    }

    // free eq_sequence 
    for (size_t i=0; i < seq_len; i++) {
        free(eq_sequence[i].coeffs);
    }
    free(eq_sequence);

    // free points
    free(points.data);

    printf("test_eq_poly_sequence_cross_check PASSED\n");
}
