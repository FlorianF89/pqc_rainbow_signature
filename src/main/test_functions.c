#include "test_functions.h"
#include <stdio.h>
#include <string.h>
#include "../lib/gf16.h"
#include "../lib/parameters.h"
#include "../lib/keygen.h"
#include "../lib/sign.h"
#include "../lib/rng.h"
#include "../lib/error_codes.h"
#include <stdlib.h>
#include <time.h>

void print_test_results(int result){
    if(result == TEST_SUCCESS){
        printf("TEST PASSED\n");
    }else{
        printf("TEST FAILED\n");
    }
}



#include "openssl/sha.h"


#ifndef _HASH_LEN
#define _HASH_LEN (32)
#endif



static inline
int _hash( unsigned char * digest , const unsigned char * m , unsigned long long mlen )
{
	SHA256_CTX sha256;
	SHA256_Init( &sha256 );
	SHA256_Update( &sha256 , m , mlen );
	SHA256_Final( digest , &sha256 );
	return 0;
}

static void print_bitsliced(const bitsliced_gf16_t in, unsigned int bit_position) {
    int is_first = 0;
    if ((in[0] >> bit_position) & 0x01u) {
        printf("1");
        is_first = 1;
    }
    if ((in[1] >> bit_position) & 0x01u) {
        if (is_first) {
            printf(" + ");
        }
        printf("y");
        is_first = 1;
    }
    if ((in[2] >> bit_position) & 0x01u) {
        if (is_first) {
            printf(" + ");
        }
        printf("X");
        is_first = 1;
    }
    if ((in[3] >> bit_position) & 0x01u) {
        if (is_first) {
            printf(" + ");
        }
        printf("y*X");
        is_first = 1;
    }
    if (is_first == 0) {
        printf("0");
    }
}


int test_bitsliced_muladd(){
    bitsliced_gf16_t i, j, i_times_j, results;
    results[0] =0;
    results[1] =0;
    results[2] =0;
    results[3] =0;
    i[0] = 0xFFFF0000FFFF0000u;
    i[1] = 0xFFFFFFFF00000000u;
    i[2] = 0;
    i[3] = 0;
    j[0] = 0xAAAAAAAAAAAAAAAAu;
    j[1] = 0xCCCCCCCCCCCCCCCCu;
    j[2] = 0xf0f0f0f0f0f0f0f0u;
    j[3] = 0xff00ff00ff00ff00u;

    i_times_j[0] = 0x6666CCCCAAAA0000u;
    i_times_j[1] = 0xAAAA6666CCCC0000u;
    i_times_j[2] = 0xff0ff00f0f00000u;
    i_times_j[3] = 0xf0f00ff0ff000000u;

    bitsliced_muladd(results, i, j);

    for (size_t i = 0; i < 4; i++)
    {
        if(i_times_j[i] != results[i]){
            return TEST_FAIL;
        }
    }
    for (size_t k = 0; k < 8800000; k++)
    {
        bitsliced_muladd(results, i, j);
    }
    
    for (size_t i = 0; i < 4; i++)
    {
        if(i_times_j[i] != results[i]){
            return TEST_FAIL;
        }
    }
    
    return TEST_SUCCESS;
}

//variable_substitution

static void print_pols(bitsliced_gf16_t pol[6], int num_of_pol){
    for (size_t i = 0; i < num_of_pol; i++)
    {
        print_bitsliced(pol[0], i);
        printf("*Z_1^2 +");
        print_bitsliced(pol[1], i);
        printf("*Z_1Z_2 +");
        print_bitsliced(pol[2], i);
        printf("*Z_1Z_3 +");
        print_bitsliced(pol[3], i);
        printf("*(Z_2)^2 +");
        print_bitsliced(pol[4], i);
        printf("*Z_2Z_3 +");
        print_bitsliced(pol[5], i);
        printf("*(Z_3)^2");
        putchar(10);
    }
}

int test_variable_substitution(){
    

    polynomial_t pol, pol1;

    bitsliced_gf16_t variable_list[O1 + V1], tmp;
    memset(variable_list, 0x00, sizeof(variable_list));
    memset(pol, 0x00, sizeof(pol));
    pol[2970][0] = 1; //f_0 = x_0^2, f_1 = ... = f_{n-1} = 0
    pol[2971][0] = 1; //f_0 = x_0^2 + x_0x_v1, f_1 = ... = f_{n-1} = 0
    variable_list[V1 + 1][0] = 0xFFFFFFFF00000000; // x0 = x0 + x_v1 + ... + v_{n-1}
    variable_list[V1][0] = 0xFFFFFFFF00000000; // x36 = x36 + x_{v1 + o1} + ... + x_{v_{n-1}}
    //for (size_t j = 0; j < 1000000; j++)
    //{
    //    polynomials_evaluation_full(tmp, pol, variable_list); 
    //    bitsliced_addition(variable_list[0], variable_list[0], tmp);
    //}
    size_t index = 0;
    for (size_t i = 0; i < NUMBER_OF_VARIABLES; i++)
    {
        int nl = 0;
        for (size_t j = i; j < NUMBER_OF_VARIABLES; j++)
        {
            if(pol[index][0] | pol[index][1] | pol[index][2] | pol[index][3]){
                printf("x_%lu x_%lu: ", i, j);
                print_bitsliced(pol[index], 0);
                putchar(' ');
                nl = 1;
            }
            index++;
        }
        if(nl)
            putchar(10);
    }
    variable_substitution_full(pol1, pol, variable_list);

    putchar(10);
    index = 0;
    for (size_t i = 0; i < NUMBER_OF_VARIABLES; i++)
    {
        int nl = 0;
        for (size_t j = i; j < NUMBER_OF_VARIABLES; j++)
        {
            if(pol1[index][0] | pol1[index][1] | pol1[index][2] | pol1[index][3]){
                printf("x_%lux_%lu: ", i, j);
                print_bitsliced(pol1[index], 0);
                putchar(' ');
                nl = 1;
            }
            index++;
        }
        if(nl)
            putchar(10);
    }
    
    return TEST_FAIL;
}

int test_solve_system(){
    
    for (size_t i = 0; i < 100; i++)
    {
        uint8_t seed[48] = {0};
        randombytes_init(seed, NULL, NULL);
        bitsliced_gf16_t sol, coef, original_sys[32], system[32], tmp;

        randombytes((uint8_t *) coef, sizeof(coef));
        randombytes((uint8_t *) system, sizeof(system));

        memcpy(original_sys, system, sizeof(system));
        solve_32x32_gf16_system(sol, system, coef);
        
        matrix_vector_multiplication(tmp, original_sys, 32, sol);
        
        if((tmp[0] & 0xFFFFFFFF) != (coef[0] & 0xFFFFFFFF) ||
            (tmp[1] & 0xFFFFFFFF) != (coef[1] & 0xFFFFFFFF) ||
            (tmp[2] & 0xFFFFFFFF) != (coef[2] & 0xFFFFFFFF) ||
            (tmp[3] & 0xFFFFFFFF) != (coef[3] & 0xFFFFFFFF)){
            return TEST_FAIL;
        }
    }
    return TEST_SUCCESS;
}

int test_signature(){
    return TEST_FAIL;
}

int test_generate_s(){
    
    bitsliced_gf16_t s[O1];
    if(generate_random_s(s) != SUCCESS){
        return TEST_FAIL;
    }
    for (size_t i = 0; i < O1; i++)
    {
        if(s[i][0] & 0xFFFFFFFF00000000 ||
           s[i][1] & 0xFFFFFFFF00000000 ||
           s[i][2] & 0xFFFFFFFF00000000 ||
           s[i][3] & 0xFFFFFFFF00000000){
               return TEST_FAIL;
        }
    }
    return TEST_SUCCESS;
}

int test_generate_t(){
    
    bitsliced_gf16_t t[O1 + V1];

    if(generate_random_t(t) != SUCCESS){
        return TEST_FAIL;
    }

    for (size_t i = 0; i < O1; i++)
    {
        if(t[i + V1][0] & 0xFFFFFFFF ||
           t[i + V1][1] & 0xFFFFFFFF ||
           t[i + V1][2] & 0xFFFFFFFF ||
           t[i + V1][3] & 0xFFFFFFFF){
               return TEST_FAIL;
        }
    }
    return TEST_SUCCESS;
}

int test_multiply_y_by_t(){
    
    bitsliced_gf16_t t[O1 + V1], variables[2], result[2], tmp;
    randombytes((uint8_t *) variables, sizeof(variables));

    variables[0][0] &= (((uint64_t) 1) << V1) - 1;
    variables[0][1] &= (((uint64_t) 1) << V1) - 1;
    variables[0][2] &= (((uint64_t) 1) << V1) - 1;
    variables[0][3] &= (((uint64_t) 1) << V1) - 1;

    //actually sets t to the identity in this representation
    memset(t, 0, sizeof(t));

    multiply_y_by_t(result, t, variables);
    bitsliced_addition(tmp, variables[0], result[0]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }
    bitsliced_addition(tmp, variables[1], result[1]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }

    //Tests with T1 T2 and T3 set to identity
    for (size_t i = 0; i < O1; i++)
    {
        t[i][0] = (((uint64_t) 1) << i) | (((uint64_t) 1) << i + O1);
    }
    for (size_t i = 0; i < O1; i++)
    {
        t[i + V1][0] = (((uint64_t) 1) << i + O1);
    }
    multiply_y_by_t(result, t, variables);
    tmp[0] = variables[0][0] ^ (variables[1][0] & 0xFFFFFFFF) ^ (variables[1][0] >> 32);
    tmp[1] = variables[0][1] ^ (variables[1][1] & 0xFFFFFFFF) ^ (variables[1][1] >> 32);
    tmp[2] = variables[0][2] ^ (variables[1][2] & 0xFFFFFFFF) ^ (variables[1][2] >> 32);
    tmp[3] = variables[0][3] ^ (variables[1][3] & 0xFFFFFFFF) ^ (variables[1][3] >> 32);

    bitsliced_addition(tmp, tmp, result[0]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }

    tmp[0] = (variables[1][0] & 0xFFFFFFFF) ^ (variables[1][0] >> 32) ^ (variables[1][0] & 0xFFFFFFFF00000000);
    tmp[1] = (variables[1][1] & 0xFFFFFFFF) ^ (variables[1][1] >> 32) ^ (variables[1][1] & 0xFFFFFFFF00000000);
    tmp[2] = (variables[1][2] & 0xFFFFFFFF) ^ (variables[1][2] >> 32) ^ (variables[1][2] & 0xFFFFFFFF00000000);
    tmp[3] = (variables[1][3] & 0xFFFFFFFF) ^ (variables[1][3] >> 32) ^ (variables[1][3] & 0xFFFFFFFF00000000);

    bitsliced_addition(tmp, tmp, result[1]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }


    return TEST_SUCCESS;
}

int test_invert_t(){
    
    bitsliced_gf16_t t[O1 + V1], t_inverse[O1 + V1], variables[2], result[2], test_values[2], tmp;

    //Tests with T1 T2 and T3 set to identity
    memset(t, 0, sizeof(t));
    for (size_t i = 0; i < O1; i++)
    {
        t[i][0] = (((uint64_t) 1) << i) | (((uint64_t) 1) << i + O1);
    }
    for (size_t i = 0; i < O1; i++)
    {
        t[i + V1][0] = (((uint64_t) 1) << i + O1);
    }

    randombytes((uint8_t *) variables, sizeof(variables));

    variables[0][0] &= (((uint64_t) 1) << V1) - 1;
    variables[0][1] &= (((uint64_t) 1) << V1) - 1;
    variables[0][2] &= (((uint64_t) 1) << V1) - 1;
    variables[0][3] &= (((uint64_t) 1) << V1) - 1;

    multiply_y_by_t(result, t, variables);
    invert_t(t_inverse, t);
    multiply_y_by_t(test_values, t_inverse, result);
    
    bitsliced_addition(tmp, variables[0], test_values[0]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }
    
    bitsliced_addition(tmp, variables[1], test_values[1]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }


    //Test with random T
    if(generate_random_t(t) != SUCCESS){
        return TEST_FAIL;
    }

    multiply_y_by_t(result, t, variables);
    invert_t(t_inverse, t);
    multiply_y_by_t(test_values, t_inverse, result);
    
    bitsliced_addition(tmp, variables[0], test_values[0]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }
    
    bitsliced_addition(tmp, variables[1], test_values[1]);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }
    
    return TEST_SUCCESS;
}

int test_variable_substitution_full(){
    
    bitsliced_gf16_t t[O1 + V1], variables[2], result[2], test_values[2], tmp, tmp1;
    polynomial_t f, f_o_t;

    randombytes((uint8_t *) variables, sizeof(variables));

    variables[0][0] &= (((uint64_t) 1) << V1) - 1;
    variables[0][1] &= (((uint64_t) 1) << V1) - 1;
    variables[0][2] &= (((uint64_t) 1) << V1) - 1;
    variables[0][3] &= (((uint64_t) 1) << V1) - 1;   
    variables[1][0] &= 0xFFFFFFFF;
    variables[1][1] &= 0xFFFFFFFF;
    variables[1][2] &= 0xFFFFFFFF;
    variables[1][3] &= 0xFFFFFFFF;    
/*     variables[0][0] = 1;
    variables[0][1] = 0;
    variables[0][2] = 0;
    variables[0][3] = 0; 
    variables[1][0] = -1lu;
    variables[1][1] = 0;
    variables[1][2] = 0;
    variables[1][3] = 0;
 */    
    //generate_random_f(f);
    memset(f, 0, sizeof(polynomial_t));
    for (size_t i = 0; i < V1; i++)
    {
        randombytes((uint8_t *)f[i], sizeof(bitsliced_gf16_t));
    }
    // randombytes((uint8_t *)f[V1 + O1], sizeof(bitsliced_gf16_t));
    generate_random_t(t);
    //memset(t, 0x00, sizeof(t));
    
    //for (size_t i = V1; i < V1 + O1; i++)
    //{
    //    t[i][0] = 0;
    //    t[i][1] = 0;
    //    t[i][2] = 0;
    //    t[i][3] = 0;
    //}

    multiply_y_by_t(test_values, t, variables);
    // put test_values and variables in right format
    shift_left_gf16(tmp, variables[1], V1);
    bitsliced_addition(variables[0], variables[0], tmp);
    shift_right_gf16(variables[1], variables[1], 64 - V1);

    shift_left_gf16(tmp, test_values[1], V1);
    bitsliced_addition(test_values[0], test_values[0], tmp);
    shift_right_gf16(test_values[1], test_values[1], 64 - V1);

    variable_substitution_full(f_o_t, f, t);

    polynomials_evaluation_full(tmp, f, test_values);
    polynomials_evaluation_full(tmp1, f_o_t, variables);
    
    bitsliced_addition(tmp, tmp, tmp1);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }

    return TEST_SUCCESS;
}