#ifndef COMPRESED_POLY_H
#define COMPRESED_POLY_H


#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../types.h"
#include "../field.h"
#include "../evaluations.h"
#include "../mle_poly.h"


void test_compress(void);
void test_sum_standard_case(void);
void test_sum_all_zero_coefficients(void);
void test_coeffs_reconstruction_standard_case(void);
void test_coeffs_all_zero_coefficients(void);
void test_coeffs_large_coefficients(void);



#endif




