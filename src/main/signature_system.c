//#include "../lib/keygen.h"
//#include "../lib/sign.h"
//#include "../lib/verify.h"
#include <stdio.h>
#include <x86intrin.h>
#include <string.h>
#include <time.h>
#include "test_functions.h"



int main(int argc, char **argv) {
    
    int total_test = 0, test_failed = 0;

    printf("Executing test_bitsliced_muladd()... ");
    int test_result = test_bitsliced_muladd();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_variable_substitution()... ");
    test_result = test_variable_substitution();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_solve_system()... ");
    test_result = test_solve_system();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_generate_s()... ");
    test_result = test_generate_s();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_generate_t()... ");
    test_result = test_generate_t();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_multiply_y_by_t()... ");
    test_result = test_multiply_y_by_t();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_invert_t()... ");
    test_result = test_invert_t();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_variable_substitution_first_layer()... ");
    test_result = test_variable_substitution_first_layer();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    printf("Executing test_variable_substitution_second_layer()... ");
    test_result = test_variable_substitution_second_layer();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);

    return 0;
}
