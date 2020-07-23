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

int generate_random_matrix_s(matrix_s_t s, prng_t *prng);

int generate_random_matrix_t(matrix_t_t t, prng_t *prng);

int generate_random_matrices_f(matrix_fi_t f[O1 + O2], prng_t *prng);

#endif //PQC_RAINBOW_KEYGEN_H
