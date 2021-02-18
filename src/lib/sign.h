//
// Created by Florian Caullery on 8/2/20.
//

#ifndef PQC_RAINBOW_SIGN_H
#define PQC_RAINBOW_SIGN_H

#include "gf16.h"
#include "utils_prng.h"
#include "keygen.h"


void print_32x32_gf16_system(bitsliced_gf16_t m[32], bitsliced_gf16_t solution);

int solve_32x32_gf16_system(bitsliced_gf16_t solution, bitsliced_gf16_t equations_coefficients[32],
                            bitsliced_gf16_t linear_coefficients);

#endif //PQC_RAINBOW_SIGN_H
