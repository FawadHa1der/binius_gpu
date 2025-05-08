#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>  // If you need NEON intrinsics
#include "unity.h"

#include "test_field.h"
#include "test_matrix_utils.h"
#include "test_mle.h"
#include "test_evaluations.h"
#include "test_compressed_poly.h"
#include "test_univariate_poly.h"
#include "test_efficient_matrix.h"
#include "test_matrix_linear.h"
#include "test_prod_check.h"
#include "test_linear_check.h"
#include "test_multi_claim.h"
#include "test_and_package.h"
#include "test_boolcheck.h"
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
    RUN_TEST(test_compute_gammas_folding_large_gamma);
    RUN_TEST(test_compute_gammas_folding_identity);
    RUN_TEST(test_compute_gammas_folding_zero_gamma);
    RUN_TEST(test_compute_gammas_folding_simple_case);

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
    RUN_TEST(test_restrict_polynomials);
    RUN_TEST(test_eq_poly);
    RUN_TEST(test_eq_poly_evaluate);
    RUN_TEST(test_eq_poly_sequence);
    RUN_TEST(test_eq_poly_sequence_random_values);
    RUN_TEST(test_evaluate_at);
    RUN_TEST(test_multilinear_lagrangian_bitand);
    RUN_TEST(test_pi);
    RUN_TEST(test_pi_all_zeroes);
    RUN_TEST(test_pi_single_non_zero);
    RUN_TEST(test_pi_alternating);
    RUN_TEST(test_pi_random_orbit);
    RUN_TEST(test_twist_untwist);
    RUN_TEST(test_twist_all_zeros);
    RUN_TEST(test_untwist_all_zeros);

    RUN_TEST(test_compress);
    RUN_TEST(test_sum_standard_case);
    RUN_TEST(test_sum_all_zero_coefficients);
    RUN_TEST(test_coeffs_reconstruction_standard_case);
    RUN_TEST(test_coeffs_all_zero_coefficients);
    RUN_TEST(test_coeffs_large_coefficients);

    // univariate polynomial tests
    RUN_TEST(test_univariate_univariate_polynomial_evaluate_at);
    RUN_TEST(test_from_evaluations_deg2);
    RUN_TEST(test_multiply_degree2_by_degree1);

    // effiecient_matrix tests    
    RUN_TEST(test_from_cols_with_queries);
    RUN_TEST(test_from_cols_single_non_zero_column);
    RUN_TEST(test_from_cols_all_zeros);
    RUN_TEST(test_from_rows_alternate_rows_full_128_bits);
    RUN_TEST(test_from_rows_all_rows_max);
    RUN_TEST(test_from_rows_single_element_in_first_row);
    RUN_TEST(test_from_rows_all_ones);
    RUN_TEST(test_from_rows_all_zeros);
    RUN_TEST(test_apply_against_traditional_matrix);
    RUN_TEST(test_frobenius_inv_lc_all_zeros);
    RUN_TEST(test_frobenius_inv_lc_single_gamma);
    RUN_TEST(test_frobenius_inv_lc_multiple_nonzero_gammas);
    RUN_TEST(test_frobenius_inv_lc_wraparound_indices);
    RUN_TEST(test_apply_all_zeros);
    RUN_TEST(test_apply_single_nonzero_byte);

    //matrix_linear tests
    RUN_TEST(test_apply_identity_matrix);
    RUN_TEST(test_apply_transposed_identity_matrix);
    RUN_TEST(test_apply_arbitrary_matrix);
    RUN_TEST(test_apply_transposed_arbitrary_matrix);
    // RUN_TEST(test_apply_invalid_input_size);
    // RUN_TEST(test_apply_invalid_output_size);


    // prod_check tests
    // RUN_TEST(test_prodcheck_new_invalid_claim); <-- SHOULD FAIL with an abort, not sure how to handle with unity since the farmework does not catch it
    RUN_TEST(test_prodcheck_new_without_checking_claim);
    RUN_TEST(test_prodcheck_new_valid_claim);
    // RUN_TEST(test_prodcheck_new_invalid_polynomial_size_p);<-- SHOULD FAIL with an abort, not sure how to handle with unity since the farmework does not catch it
    // RUN_TEST(test_prodcheck_new_invalid_polynomial_size_q); //<-- SHOULD FAIL with an abort, not sure how to handle with unity since the farmework does not catch it
    RUN_TEST(test_prodcheck_compute_round_polynomial_valid);
    //RUN_TEST(test_prodcheck_compute_round_polynomial_invalid_claim); // SHOULD FAIL
    // RUN_TEST(test_prodcheck_compute_round_polynomial_protocol_complete); // SHOULD FAIL
    RUN_TEST(test_prodcheck_bind_valid_challenge);
    RUN_TEST(test_prodcheck_finish_valid_state);
    // RUN_TEST(test_prodcheck_finish_incomplete_state);// SHOULD FAIL
    RUN_TEST(test_prodcheck_finish_multiple_variables);
    RUN_TEST(test_prodcheck_full);


    // linear_check tests
    RUN_TEST(test_new_lincheck);
    //RUN_TEST(test_invalid_matrix_dimensions); // SHOULD FAIL
    //RUN_TEST(test_invalid_polynomial_length); // SHOULD FAIL
    RUN_TEST(test_lincheck_builder_new_with_valid_inputs);
    RUN_TEST(test_lincheck_builder_build);
    //RUN_TEST(test_lincheck_builder_invalid_matrix_input_size); // SHOULD FAIL
    //RUN_TEST(test_lincheck_builder_invalid_matrix_output_size); // SHOULD FAIL
    //RUN_TEST(test_lincheck_builder_insufficient_points); // SHOULD FAIL
    //RUN_TEST(test_lincheck_builder_invalid_polynomial_length); // SHOULD FAIL

    // multi_claim tests
    //RUN_TEST(test_multiclaim_builder_default);
    RUN_TEST(test_multiclaim_builder_new_valid);
    //RUN_TEST(test_multiclaim_builder_new_invalid_polynomial_length); // SHOULD FAIL
    RUN_TEST(test_multiclaim_builder_new_edge_case);
    RUN_TEST(test_multiclaim_builder_build_simple_case);
    RUN_TEST(test_new_default_inputs);
    RUN_TEST(test_new_with_nonzero_openings);
    // RUN_TEST(test_multiclaim_complete); // TAKE A LONG TIME, TODO UNCOMMENT WHEN READY


    // and_package tests
    RUN_TEST(test_exec_alg_and);

    //boolcheck tests
    RUN_TEST(test_trit_mapping_small_c);
    RUN_TEST(test_trit_mapping_medium_c);
    RUN_TEST(test_trit_mapping_large_c);
    RUN_TEST(test_trit_mapping_no_c);
    RUN_TEST(test_extend_n_tables);
    



    return UNITY_END();
}
