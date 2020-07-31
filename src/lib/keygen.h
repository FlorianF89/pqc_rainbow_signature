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

typedef struct private_key {
    matrix_s_t s;
    matrix_t_t t;
    matrix_fi_t f[O1 + O2];
    gf16_32x32_matrix s_prime;
    gf16_32x32_matrix t1;
    gf16_32x32_matrix t2;
    gf16_32x32_matrix t3;
    gf16_32x32_matrix f1s[N - V1];
    gf16_32x32_matrix f2s[N - V1];
    gf16_32x32_matrix f3s[N - V2];
    gf16_32x32_matrix f5s[N - V2];
    gf16_32x32_matrix f6s[N - V2];
} private_key_t;

int generate_private_key(private_key_t *private_key, prng_t *prng);

int generate_random_matrix_s(private_key_t *private_key, prng_t *prng);

int generate_random_matrix_t(private_key_t *private_key, prng_t *prng);

int generate_random_matrices_f(private_key_t *private_key, prng_t *prng);

void transpose_and_add_32x32_gf16_matrices(matrix_s_t s1_plus_s2t, matrix_s_t s1, matrix_s_t s2);

void gaussian_elimination_for_32x32_gf16_matrix(matrix_s_t s);

void set_32x32_gf16_matrix_to_identity(bitsliced_gf16_t a[32]);

void multiply_32x32_gf16_matrices(bitsliced_gf16_t a_times_b[32], bitsliced_gf16_t a[32], bitsliced_gf16_t b[32]);

#endif //PQC_RAINBOW_KEYGEN_H
