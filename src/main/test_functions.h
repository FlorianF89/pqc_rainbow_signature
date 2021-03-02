#ifndef MAIN_HELLO_GREET_H_
#define MAIN_HELLO_GREET_H_

#define TEST_SUCCESS 0
#define TEST_FAIL (-1)

int test_bitsliced_muladd();

int test_solve_system();

int test_generate_s();

int test_generate_t();

int test_multiply_y_by_t();

int test_invert_t();

int test_variable_substitution_first_layer();

int test_variable_substitution_second_layer();

int test_multiply_by_s();

int test_variable_private_key_to_public_key();

int test_find_preimage_by_central_map();

int test_signature();

int test_verify_positive();

#endif
