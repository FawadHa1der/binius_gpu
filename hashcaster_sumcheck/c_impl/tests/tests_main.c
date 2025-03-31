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

    RUN_TEST(test_drop_top_bit_standard_cases);
    RUN_TEST(test_drop_top_bit_edge_cases);
    RUN_TEST(test_drop_top_bit_large_numbers);
    RUN_TEST(test_drop_top_bit_all_bits_set);

    RUN_TEST(test_cpu_v_movemask_epi8_standard_cases);
    RUN_TEST(test_cpu_v_movemask_epi8_all_ones);
    RUN_TEST(test_cpu_v_movemask_epi8_all_zeros);

    RUN_TEST(test_v_slli_epi64_basic_shift);
    RUN_TEST(test_v_slli_epi64_zero_shift);
    RUN_TEST(test_v_slli_epi64_edge_cases);
    RUN_TEST(test_v_slli_epi64_no_overflow);



    return UNITY_END();
}
