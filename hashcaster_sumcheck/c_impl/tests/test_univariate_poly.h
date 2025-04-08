#ifndef UNIVARIATE_POLY_TEST_H
#define UNIVARIATE_POLY_TEST_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include "../field.h"
#include "../types.h"
#include "../univariate_poly.h"


void test_univariate_polynomial_evaluate_at(void);
void test_from_evaluations_deg2(void);
void test_multiply_degree2_by_degree1(void);
#endif