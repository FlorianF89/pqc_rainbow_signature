//
// Created by Florian Caullery on 7/20/20.
//

#ifndef PQC_RAINBOW_KEYGEN_H
#define PQC_RAINBOW_KEYGEN_H

#include "gf16.h"
#include "parameters.h"
#include "utils_prng.h"

#define ELMTS_IN_BITSLICED 64

typedef bitsliced_gf16_t first_layer_polynomial_t[V1][2];

typedef bitsliced_gf16_t second_layer_polynomial_t[V1 + O1][2];

typedef struct equation_system_t{
    first_layer_polynomial_t first_layer_polynomials[O1];
    second_layer_polynomial_t second_layer_polynomials[O2];
}equation_system_t;

void variable_substitution_first_layer(bitsliced_gf16_t f_times_t[NUMBER_OF_VARIABLES][2], 
                                            first_layer_polynomial_t original, 
                                            bitsliced_gf16_t variable_list[V1 + O1]);

void variable_substitution_second_layer(bitsliced_gf16_t f_times_t[NUMBER_OF_VARIABLES][2], 
                                            second_layer_polynomial_t original, 
                                            bitsliced_gf16_t variable_list[V1 + O1]);

int variable_substitution_full(equation_system_t transformed, equation_system_t original, bitsliced_gf16_t variable_list[V1 + O1]);

void polynomial_evaluation_first_layer(bitsliced_gf16_t result, first_layer_polynomial_t pol, 
                                bitsliced_gf16_t variables_value[2], int position);

void polynomial_evaluation_second_layer(bitsliced_gf16_t results, first_layer_polynomial_t pol, 
                                bitsliced_gf16_t variables_value[2], int position);

void polynomial_evaluation_full(bitsliced_gf16_t results, bitsliced_gf16_t pol[NUMBER_OF_VARIABLES][2], 
                                bitsliced_gf16_t variables_value[2], int position);

void scalar_product(uint8_t pos_to_put_res, bitsliced_gf16_t res, bitsliced_gf16_t in1, bitsliced_gf16_t in2);

void matrix_vector_multiplication(bitsliced_gf16_t res, bitsliced_gf16_t *mat, int matrix_width, bitsliced_gf16_t vec);

void matrix_matrix_multiplication(bitsliced_gf16_t *res, 
                                bitsliced_gf16_t *mat_a, int a_width, 
                                bitsliced_gf16_t *mat_b, int b_width);

int generate_random_s(bitsliced_gf16_t s[O1]);

int generate_random_t(bitsliced_gf16_t t[O1 + V1]);

void multiply_y_by_t(bitsliced_gf16_t result[2], bitsliced_gf16_t t[V1 + O1],
                     bitsliced_gf16_t variables[2]);

void invert_t(bitsliced_gf16_t t_inverse[V1 + O1], bitsliced_gf16_t t[V1 + O1]);

//int generate_random_f(polynomial_t f);

#endif //PQC_RAINBOW_KEYGEN_H
