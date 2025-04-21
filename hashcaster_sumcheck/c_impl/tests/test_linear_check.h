#ifndef TEST_LINEARCHECK_H
#define TEST_LINEARCHECK_H


void test_default_lincheck(void);
void test_new_lincheck(void);
void test_invalid_matrix_dimensions();
void test_invalid_polynomial_length();
void test_lincheck_builder_new_with_valid_inputs();
void test_lincheck_builder_build(void);
void test_lincheck_builder_invalid_matrix_input_size(void);
void test_lincheck_builder_invalid_matrix_output_size(void);
void test_lincheck_builder_insufficient_points(void) ;
void test_lincheck_builder_invalid_polynomial_length(void) ;
#endif