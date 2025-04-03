#ifndef TEST_MLE_H
#define TEST_MLE_H


#include "../mle_poly.h"

void test_eq_eval_basic(void);
void test_eq_eval_identity(void);
void test_eq_eval_empty(void);
void test_eq_poly_sequence_cross_check(void);
void test_to_points_inv_orbit_ones(void);
void test_to_points_inv_orbit_last_element(void);
void test_eq_sums(void) ;

void test_drop_top_bit_standard_cases(void);
void test_drop_top_bit_zero(void);

void test_drop_top_bit_edge_cases(void);
void test_drop_top_bit_large_numbers(void);
void test_drop_top_bit_all_bits_set(void);




void test_cpu_v_movemask_epi8_standard_cases(void);

void test_cpu_v_movemask_epi8_all_ones(void);
void test_cpu_v_movemask_epi8_all_zeros(void);
void test_v_slli_epi64_basic_shift(void);
void test_v_slli_epi64_zero_shift(void);
void test_v_slli_epi64_edge_cases(void);
void test_v_slli_epi64_no_overflow(void);
void test_restrict_polynomials(void);
void test_eq_poly(void);
void test_eq_poly_evaluate(void);
void test_eq_poly_sequence(void);
void test_eq_poly_sequence_random_values(void);
void test_evaluate_at(void);
void test_multilinear_lagrangian_bitand(void);
#endif