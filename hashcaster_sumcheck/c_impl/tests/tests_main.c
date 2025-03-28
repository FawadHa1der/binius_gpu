#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>  // If you need NEON intrinsics
#include "unity.h"

#include "test_field.h"
#include "test_matrix_utils.h"
#include "test_mle.h"
// Unity expects these, even if empty
void setUp(void)
{
    // optional: do something before each test
}

void tearDown(void)
{
    // optional: do something after each test
}

/*-------------------------------------------------------
 * Unity's main test runner
 --------------------------------------------------------*/
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_f128_is_field);
    RUN_TEST(test_frobenius);
    RUN_TEST(test_pi_as_expected);
    RUN_TEST(test_twists_logic_and);

    RUN_TEST(test_apply_matrix);
    RUN_TEST(test_invert_matrix);

    // MLE tests
    RUN_TEST(test_eq_eval_basic);
    RUN_TEST(test_eq_eval_identity);
    RUN_TEST(test_eq_eval_empty);
    RUN_TEST(test_eq_poly_sequence_cross_check);
    RUN_TEST(test_to_points_inv_orbit_ones);
    RUN_TEST(test_to_points_inv_orbit_last_element);
    RUN_TEST(test_eq_sums);

    return UNITY_END();
}
