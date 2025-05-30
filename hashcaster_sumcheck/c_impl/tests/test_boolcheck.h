#ifndef TEST_BOOLCHECK_H
#define TEST_BOOLCHECK_H


#include "../boolcheck.h"
#include "../field.h"


void test_trit_mapping_small_c(void) ;
void test_trit_mapping_medium_c(void) ;
void test_trit_mapping_large_c(void) ;
void test_trit_mapping_no_c(void) ;
void test_extend_n_tables(void) ;
void test_new_andcheck1(void) ;
void test_new_andcheck_with_multiclaim(void) ;
void test_compute_imaginary_rounds(void) ;
void test_compute_round_polynomial_cached_result(void) ;
void test_compute_round_polynomial_initial_round(void) ;
void test_compute_round_polynomial_exceeds_round_limit(void) ;
void test_compute_round_polynomial_invalid_claim(void) ;
void test_compute_round_polynomial_correct_claim_update(void) ;
void test_bind_single(void) ;


#endif