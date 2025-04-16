#ifndef TEST_PROD_CHECK_H
#define TEST_PROD_CHECK_H


#include "../prod_check.h"
#include "../field.h"
#include "../mle_poly.h"
#include "../compressed_poly.h"



void test_prodcheck_new_invalid_claim(void);
void test_prodcheck_new_without_checking_claim(void);
void test_prodcheck_new_valid_claim(void);
void test_prodcheck_new_invalid_polynomial_size_p(void);
void test_prodcheck_new_invalid_polynomial_size_q(void);
void test_prodcheck_compute_round_polynomial_valid(void);
void test_prodcheck_compute_round_polynomial_invalid_claim(void);
void test_prodcheck_compute_round_polynomial_protocol_complete(void);
void test_prodcheck_bind_valid_challenge(void);
void test_prodcheck_finish_valid_state(void);
void test_prodcheck_finish_incomplete_state(void);
void test_prodcheck_finish_multiple_variables(void);
void test_prodcheck_full(void);

#endif