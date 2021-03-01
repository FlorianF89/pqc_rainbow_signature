//#include "../lib/keygen.h"
//#include "../lib/sign.h"
//#include "../lib/verify.h"
#include <stdio.h>
#include <x86intrin.h>
#include <string.h>
#include <time.h>
#include "test_functions.h"


#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define RESET "\x1B[0m"

void print_test_results(int result){
    if(result == TEST_SUCCESS){
        printf(GRN "PASSED:\t");
    }else{
        printf(RED "FAILED:\t");
    }
}

int main(int argc, char **argv) {
    
    int total_test = 0, test_failed = 0;

    int test_result = test_bitsliced_muladd();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_bitsliced_muladd()\n" RESET);

    test_result = test_solve_system();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_solve_system()\n" RESET);

    test_result = test_generate_s();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_generate_s()\n" RESET);

    test_result = test_generate_t();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_generate_t()\n" RESET);

    test_result = test_multiply_y_by_t();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_multiply_y_by_t()\n" RESET);

    test_result = test_invert_t();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_invert_t()\n" RESET);

    test_result = test_variable_substitution_first_layer();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_variable_substitution_first_layer()\n" RESET);

    test_result = test_variable_substitution_second_layer();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_variable_substitution_second_layer()\n" RESET);

    test_result = test_multiply_by_s();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_multiply_by_s()\n" RESET);

    test_result = test_variable_private_key_to_public_key();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_variable_private_key_to_public_key()\n" RESET);

    test_result = test_find_preimage_by_central_map();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_find_preimage_by_central_map()\n" RESET);

    test_result = test_signature();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_signature()\n" RESET);

    test_result = test_verify_positive();
    total_test++;
    test_failed += test_result == TEST_SUCCESS ? 0 : 1;
    print_test_results(test_result);
    printf("test_verify()\n" RESET);
    return 0;
}
