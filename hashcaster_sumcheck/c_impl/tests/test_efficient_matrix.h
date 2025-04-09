#ifndef TEST_EFFICIENT_MATRIX_H
#define TEST_EFFICIENT_MATRIX_H
#include "../efficient_matrix.h"
#include "../types.h"
#include "../field.h"


void test_from_cols_single_non_zero_column(void);
void test_from_cols_all_zeros(void);
void test_from_cols_with_queries(void);
void test_from_rows_all_zeros(void);
void test_from_rows_all_ones(void);
void test_from_rows_single_element_in_first_row(void);
void test_from_rows_all_rows_max(void);
void test_from_rows_alternate_rows_full_128_bits(void);
void test_apply_all_zeros(void);
void test_apply_single_nonzero_byte(void);
void test_apply_against_traditional_matrix(void);
void test_frobenius_inv_lc_all_zeros(void);
void test_frobenius_inv_lc_single_gamma(void);
void test_frobenius_inv_lc_multiple_nonzero_gammas(void);
void test_frobenius_inv_lc_wraparound_indices(void);

#endif