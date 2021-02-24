#ifndef MAIN_HELLO_GREET_H_
#define MAIN_HELLO_GREET_H_

#define TEST_SUCCESS 0
#define TEST_FAIL (-1)

void print_test_results(int result);

int test_bitsliced_muladd();

int test_variable_substitution();

int test_solve_system();

int test_generate_s();

int test_generate_t();

int test_multiply_y_by_t();

int test_invert_t();

int test_variable_substitution_first_layer();

int test_variable_substitution_second_layer();

#endif
