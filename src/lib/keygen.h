//
// Created by Florian Caullery on 7/20/20.
//

#ifndef PQC_RAINBOW_KEYGEN_H
#define PQC_RAINBOW_KEYGEN_H

#include "gf16.h"
#include "parameters.h"
#include "utils_prng.h"

typedef bitsliced_gf16_t matrix_s_t[O1 + O2];
typedef bitsliced_gf16_t matrix_t_t[(V1 + O1 + O2 + sizeof(uint64_t) - 1) / sizeof(uint64_t)][V1 + O1 + O2];
typedef matrix_t_t matrix_fi_t;
typedef bitsliced_gf16_t gf16_32x32_matrix[32];

typedef struct bitsliced_quadratic_polynomials {
    bitsliced_gf16_t coefficients[(N * (N + 1) / 2)];
} bitsliced_quadratic_polynomials_t;

typedef struct private_key {
    matrix_s_t s;
    matrix_t_t t;
    matrix_fi_t f[O1 + O2];
    gf16_32x32_matrix s_prime;
    bitsliced_gf16_t t1_t2[N - V1];
    bitsliced_gf16_t new_t3[N - V2];
    bitsliced_gf16_t t1_t3_minus_t2[N - V2];
    bitsliced_quadratic_polynomials_t mq_polynomials;
    gf16_32x32_matrix t1;
    gf16_32x32_matrix t2;
    gf16_32x32_matrix t3;
    gf16_32x32_matrix f1s[N - V1];
    gf16_32x32_matrix f2s[N - V1];
    gf16_32x32_matrix f3s[N - V2];
    gf16_32x32_matrix f5s[N - V2];
    gf16_32x32_matrix f6s[N - V2];
} private_key_t;


typedef struct public_key {
    bitsliced_quadratic_polynomials_t polynomials;
    bitsliced_gf16_t mq[(N * (N + 1) / 2)];
    bitsliced_gf16_t mp[(N * (N + 1) / 2)];
} public_key_t;

int generate_private_key(private_key_t *private_key, prng_t *prng);

int generate_random_matrix_s(private_key_t *private_key, prng_t *prng);

int generate_random_matrix_t(private_key_t *private_key, prng_t *prng);

int generate_random_matrices_f(private_key_t *private_key, prng_t *prng);

int generate_random_32x32_gf16_upper_triangular_matrix(bitsliced_gf16_t a[32], prng_t *prng);

int generate_random_32x32_gf16_matrix(bitsliced_gf16_t a[32], prng_t *prng);

void transpose_32x32_gf16_matrix(bitsliced_gf16_t out[32], bitsliced_gf16_t in[32]);

void bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(bitsliced_gf16_t *out,
                                                                             bitsliced_gf16_t *in,
                                                                             uint32_t position_to_place);

void bitsliced_gf16_sum_64_first_elements_and_place_result_in_given_position(bitsliced_gf16_t *out,
                                                                             bitsliced_gf16_t *in,
                                                                             uint32_t position_to_place);

uint32_t parity_of_32_bit_words(uint32_t in);

void transpose_and_add_32x32_gf16_matrices(matrix_s_t s1_plus_s2t, matrix_s_t s1, matrix_s_t s2);

void gaussian_elimination_for_32x32_gf16_matrix(matrix_s_t s);

bitsliced_gf16_t extract_then_expand(bitsliced_gf16_t a, unsigned int position, unsigned int shift);

void set_32x32_gf16_matrix_to_identity(bitsliced_gf16_t a[32]);

void multiply_32x32_gf16_matrices(bitsliced_gf16_t a_times_b[32], bitsliced_gf16_t a[32], bitsliced_gf16_t b[32]);

void derive_public_key_from_private_key(public_key_t *public_key, private_key_t *private_key);


void replace_variable_by_linear_combination_in_quadratic_polynomial(bitsliced_quadratic_polynomials_t *modified_f,
                                                                    bitsliced_quadratic_polynomials_t *original_f,
                                                                    int variable_index,
                                                                    bitsliced_gf16_t *linear_combination);

#endif //PQC_RAINBOW_KEYGEN_H
