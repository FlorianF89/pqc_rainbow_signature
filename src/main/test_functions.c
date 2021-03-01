#include "test_functions.h"
#include <stdio.h>
#include <string.h>
#include "../lib/gf16.h"
#include "../lib/parameters.h"
#include "../lib/keygen.h"
#include "../lib/sign.h"
#include "../lib/rng.h"
#include "../lib/error_codes.h"
#include "../lib/utils_hash.h"
#include <stdlib.h>
#include <time.h>



#include "openssl/sha.h"


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

int test_solve_system(){
    
    for (size_t i = 0; i < 1; i++)
    {
        uint8_t seed[48] = {0};
        randombytes_init(seed, NULL, NULL);
        bitsliced_gf16_t sol, coef, original_sys[32], system[32], tmp;

        randombytes((uint8_t *) coef, sizeof(coef));
        randombytes((uint8_t *) system, sizeof(system));
        for (size_t i = 0; i < 32; i++)
        {
            shift_right_gf16(system[i], system[i], 32);
        }

        for (int i = 0; i < 32; i++)
        {
            system[i][0] ^= ((coef[0] >> i) & 0x01llu) << 32;
            system[i][1] ^= ((coef[1] >> i) & 0x01llu) << 32;
            system[i][2] ^= ((coef[2] >> i) & 0x01llu) << 32;
            system[i][3] ^= ((coef[3] >> i) & 0x01llu) << 32;
        }

        memcpy(original_sys, system, sizeof(system));
        solve_32x32_gf16_system(sol, system);
        
        line_matrix_vector_multiplication(tmp, original_sys, 32, sol);
        
        if((tmp[0] & 0xFFFFFFFF) != (coef[0] & 0xFFFFFFFF) ||
            (tmp[1] & 0xFFFFFFFF) != (coef[1] & 0xFFFFFFFF) ||
            (tmp[2] & 0xFFFFFFFF) != (coef[2] & 0xFFFFFFFF) ||
            (tmp[3] & 0xFFFFFFFF) != (coef[3] & 0xFFFFFFFF)){
            return TEST_FAIL;
        }
    }
    return TEST_SUCCESS;
}

int test_find_preimage_by_central_map(){

    central_map_t f;
    generate_random_central_map(&f);

    bitsliced_gf16_t image, preimage[2], tmp;
    randombytes((uint8_t *) image, sizeof(bitsliced_gf16_t));

    find_preimage_by_central_map(preimage, f, image);
    evaluate_central_map(tmp, &f, preimage);

    bitsliced_addition(tmp, tmp, image);
    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }
    else{
        return TEST_SUCCESS;
    }
}

int test_generate_s(){
    
    bitsliced_gf16_t s[O1];
    if(generate_random_s(s) != SUCCESS){
        return TEST_FAIL;
    }
    for (size_t i = 0; i < O1; i++)
    {
        if(s[i][0] & 0xFFFFFFFF ||
           s[i][1] & 0xFFFFFFFF ||
           s[i][2] & 0xFFFFFFFF ||
           s[i][3] & 0xFFFFFFFF){
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

void print_first_layer_pol(first_layer_polynomial_t pol){

    for (size_t i = 0; i < V1; i++)
    {    
        printf("[");
        for (size_t j = 0; j < 36; j++)
        {
            print_bitsliced(pol[i][0], j);
            printf(", ");
        }
        for (size_t j = 0; j < 63; j++)
        {
            print_bitsliced(pol[i][1], j);
            printf(", ");
        }
        print_bitsliced(pol[i][1], 63);
        printf("],");
        putchar(10);
    }
    putchar(10);
}

void print_second_layer_pol(second_layer_polynomial_t pol){

    for (size_t i = 0; i < V1 + O1; i++)
    {    
        printf("[");
        for (size_t j = 0; j < V1; j++)
        {
            print_bitsliced(pol[i][0], j);
            printf(", ");
        }
        for (size_t j = 0; j < O1 + O2 - 1; j++)
        {
            print_bitsliced(pol[i][1], j);
            printf(", ");
        }
        print_bitsliced(pol[i][1], O1 + O2 - 1);
        printf("],");
        putchar(10);
    }
    putchar(10);
}

void print_t(bitsliced_gf16_t t[V1 + O1]){

    for (size_t i = 0; i < V1; i++)
    {    
        printf("[");
        for (size_t j = 0; j < V1; j++)
        {
            if(i == j){
                printf("1, ");
            }
            else{
                printf("0, ");
            }
        }
        for (size_t j = 0; j < 63; j++)
        {
            print_bitsliced(t[i], j);
            printf(", ");
        }
        print_bitsliced(t[i], 63);
        printf("],");
        putchar(10);
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {    
        printf("[");
        for (size_t j = 0; j < V1 + O1; j++)
        {
            if(i == j){
                printf("1, ");
            }
            else{
                printf("0, ");
            }
        }
        for (size_t j = V1 + O1; j < V1 + O1 + O2 - 1; j++)
        {
            print_bitsliced(t[i], j - V1);
            printf(", ");
        }
        print_bitsliced(t[i], 63);
        printf("],");
        putchar(10);
    }
    for (size_t i = V1 + O1; i < V1 + O1 + O2; i++)
    {    
        printf("[");
        for (size_t j = 0; j < V1 + O1 + O2 - 1; j++)
        {
            if(i == j){
                printf("1, ");
            }
            else{
                printf("0, ");
            }
        }
        if(i == 99){
                printf("1],");
            }
            else{
                printf("0],");
            }
        putchar(10);
    }
    putchar(10);
}

int test_variable_substitution_first_layer(){
    
    first_layer_polynomial_t f;
    bitsliced_gf16_t transformed[NUMBER_OF_VARIABLES][2];

    for (size_t i = 0; i < V1; i++)
    {
        randombytes((uint8_t *) &f[i][0][0], (V1 + 7) / 8);
        randombytes((uint8_t *) &f[i][0][1], (V1 + 7) / 8);
        randombytes((uint8_t *) &f[i][0][2], (V1 + 7) / 8);
        randombytes((uint8_t *) &f[i][0][3], (V1 + 7) / 8);
        f[i][0][0] = (f[i][0][0] << i) & (1llu << (V1)) - 1;
        f[i][0][1] = (f[i][0][1] << i) & (1llu << (V1)) - 1;
        f[i][0][2] = (f[i][0][2] << i) & (1llu << (V1)) - 1;
        f[i][0][3] = (f[i][0][3] << i) & (1llu << (V1)) - 1;
        randombytes((uint8_t *) &f[i][1][0], (O1 + 7) / 8);
        randombytes((uint8_t *) &f[i][1][1], (O1 + 7) / 8);
        randombytes((uint8_t *) &f[i][1][2], (O1 + 7) / 8);
        randombytes((uint8_t *) &f[i][1][3], (O1 + 7) / 8);
        f[i][1][0] &= (1llu << O1) - 1;
        f[i][1][1] &= (1llu << O1) - 1;
        f[i][1][2] &= (1llu << O1) - 1;
        f[i][1][3] &= (1llu << O1) - 1;
    }

    bitsliced_gf16_t t[O1 + V1], tmp, tmp1, variables[2], var_times_t[2];
    generate_random_t(t);    
 
    randombytes((uint8_t *) &variables[0][0], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][1], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][2], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][3], (V1 + 7) / 8);

    variables[0][0] &= (1llu << (V1)) - 1;
    variables[0][1] &= (1llu << (V1)) - 1;
    variables[0][2] &= (1llu << (V1)) - 1;
    variables[0][3] &= (1llu << (V1)) - 1;
    randombytes((uint8_t *) variables[1], sizeof(bitsliced_gf16_t));


    variable_substitution_first_layer(transformed, f, t);
    polynomial_evaluation_full(tmp, transformed, variables, 0);

    multiply_y_by_t(var_times_t, t, variables);
    polynomial_evaluation_first_layer(tmp1, f, var_times_t, 0);
    bitsliced_addition(tmp, tmp, tmp1);
    if(gf16_is_zero(tmp, 0)){
        return TEST_SUCCESS;
    }
    else{
        return TEST_FAIL;
    }
    
}

int test_variable_substitution_second_layer(){
    
    second_layer_polynomial_t f;
    bitsliced_gf16_t transformed[NUMBER_OF_VARIABLES][2];
    memset(f, 0, sizeof(f));
    for (size_t i = 0; i < V1; i++)
    {
        randombytes((uint8_t *) &f[i][0][0], (V1 + 7) / 8);
        randombytes((uint8_t *) &f[i][0][1], (V1 + 7) / 8);
        randombytes((uint8_t *) &f[i][0][2], (V1 + 7) / 8);
        randombytes((uint8_t *) &f[i][0][3], (V1 + 7) / 8);
        randombytes((uint8_t *) f[i][1], sizeof(bitsliced_gf16_t));
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {
        randombytes((uint8_t *) f[i][1], sizeof(bitsliced_gf16_t));        
    }

    bitsliced_gf16_t t[O1 + V1], tmp, tmp1, variables[2], var_times_t[2];
    generate_random_t(t);   
    memset(t, 0, sizeof(t));
    for (size_t i = 0; i < V1; i++)
    {
        t[i][0] = 0x100000001lu << i;
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {
        t[i][0] = 0x100000000lu << (i - V1);
    }

    randombytes((uint8_t *) &variables[0][0], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][1], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][2], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][3], (V1 + 7) / 8);

    variables[0][0] &= (1llu << (V1)) - 1;
    variables[0][1] &= (1llu << (V1)) - 1;
    variables[0][2] &= (1llu << (V1)) - 1;
    variables[0][3] &= (1llu << (V1)) - 1;
    randombytes((uint8_t *) variables[1], sizeof(bitsliced_gf16_t));
    variable_substitution_second_layer(transformed, f, t);
    polynomial_evaluation_full(tmp, transformed, variables, 0);

    multiply_y_by_t(var_times_t, t, variables);
    polynomial_evaluation_second_layer(tmp1, f, var_times_t, 0);
    bitsliced_addition(tmp, tmp, tmp1);
    if(gf16_is_zero(tmp, 0)){
        return TEST_SUCCESS;
    }
    else{
        return TEST_FAIL;
    }
}

int test_multiply_by_s(){
    
    //S is an involution so we test multiplying twice by S and check if the entry is unmodified
    bitsliced_gf16_t original, result, tmp, s[O1];
    randombytes((uint8_t *) original, sizeof(bitsliced_gf16_t));

    randombytes((uint8_t *) s, sizeof(s));
    for (size_t i = 0; i < O1; i++)
    {
        shift_left_gf16(s[i], s[i], O1);
    }
    
    multiply_y_by_s(tmp, s, original);
    multiply_y_by_s(result, s, tmp);

    bitsliced_addition(tmp, result, original);

    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }else{
        return TEST_SUCCESS;
    }
}

int test_variable_private_key_to_public_key(){
    central_map_t f;
    bitsliced_gf16_t public_key[O1 + O2][NUMBER_OF_VARIABLES][2], s[O1], t[O1 + V1], 
                tmp, tmp1, tmp2, variables[2], var_times_t[2];;

    generate_random_central_map(&f);    
    generate_random_t(t);    
    generate_random_s(s);
    randombytes((uint8_t *) &variables[0][0], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][1], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][2], (V1 + 7) / 8);
    randombytes((uint8_t *) &variables[0][3], (V1 + 7) / 8);

    variables[0][0] &= (1llu << (V1)) - 1;
    variables[0][1] &= (1llu << (V1)) - 1;
    variables[0][2] &= (1llu << (V1)) - 1;
    variables[0][3] &= (1llu << (V1)) - 1;
    randombytes((uint8_t *) variables[1], sizeof(bitsliced_gf16_t));


    private_key_to_public_key(public_key, f, t, s);

    memset(tmp1, 0, sizeof(bitsliced_gf16_t));
    for (size_t i = 0; i < O1 + O2; i++)
    {
        polynomial_evaluation_full(tmp, public_key[i], variables, i);
        bitsliced_addition(tmp1, tmp, tmp1);
    }

    memset(tmp2, 0, sizeof(bitsliced_gf16_t));
    multiply_y_by_t(var_times_t, t, variables);
    evaluate_central_map(tmp2, &f, var_times_t);

    multiply_y_by_s(tmp, s, tmp2);
    bitsliced_addition(tmp, tmp, tmp1);

    if(tmp[0] || tmp[1] || tmp[2] || tmp[3]){
        return TEST_FAIL;
    }else{
        return TEST_SUCCESS;
    }
}


int test_signature(){

    central_map_t f;
    generate_random_central_map(&f);

    bitsliced_gf16_t t[O1 + V1], t_inverse[O1 + V1], s[O1], tmp, tmp1[2], signature[2], public_key[O1 + O2][NUMBER_OF_VARIABLES][2];
    generate_random_t(t);
    generate_random_s(s);
    invert_t(t_inverse, t);

    uint8_t message[32] = {0};
    sign_message(signature, &f, t_inverse, s, message, sizeof(message));

    multiply_y_by_t(tmp1, t, signature);
    evaluate_central_map(tmp, &f, tmp1);
    multiply_y_by_s(tmp1[0], s, tmp);

    uint8_t digest[_HASH_LEN];
    hash_msg(digest, _HASH_LEN, message, sizeof(message));
    memcpy(tmp, digest, sizeof(bitsliced_gf16_t));

    if(tmp[0] != tmp1[0][0] || 
        tmp[1] != tmp1[0][1] || 
        tmp[2] != tmp1[0][2] || 
        tmp[3] != tmp1[0][3]){
        return TEST_FAIL;
    }
    else{
        return TEST_SUCCESS;
    }
}


int test_verify_positive(){

    central_map_t f;
    generate_random_central_map(&f);

    bitsliced_gf16_t t[O1 + V1], t_inverse[O1 + V1], s[O1], tmp, tmp1[2], signature[2], public_key[O1 + O2][NUMBER_OF_VARIABLES][2];
    generate_random_t(t);
    generate_random_s(s);
    invert_t(t_inverse, t);

    uint8_t message[32] = {0};
    sign_message(signature, &f, t_inverse, s, message, sizeof(message));
    veri

    if(tmp[0] != tmp1[0][0] || 
        tmp[1] != tmp1[0][1] || 
        tmp[2] != tmp1[0][2] || 
        tmp[3] != tmp1[0][3]){
        return TEST_FAIL;
    }
    else{
        return TEST_SUCCESS;
    }
}
