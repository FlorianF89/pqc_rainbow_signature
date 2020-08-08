//
// Created by Florian Caullery on 8/2/20.
//

#ifndef PQC_RAINBOW_SIGN_H
#define PQC_RAINBOW_SIGN_H

#include "gf16.h"
#include "utils_prng.h"

void evaluate_quadratic_polynomial_at_x0_x31(bitsliced_gf16_t *evaluation, uint8_t position_to_put_evaluation,
                                             bitsliced_gf16_t *f, bitsliced_gf16_t *x0_x31);

void evaluate_32_quadratic_polynomials_at_x0_x31(bitsliced_gf16_t *evaluation, bitsliced_gf16_t f[32][32],
                                                 bitsliced_gf16_t *x0_x31);


void evaluate_32_quadratic_polynomials_at_x0_x63(bitsliced_gf16_t *evaluation, bitsliced_gf16_t f[32][64],
                                                 bitsliced_gf16_t *x0_x63);

int solve_32x32_gf16_system(bitsliced_gf16_t *solution, bitsliced_gf16_t equations_coefficients[32],
                            bitsliced_gf16_t *linear_coefficients);

void find_preimage_of_x0_x31_by_32_polynomials_in_64_variables(bitsliced_gf16_t *preimage, bitsliced_gf16_t f[32][64],
                                                               bitsliced_gf16_t *x0_x31, prng_t *prng);


#endif //PQC_RAINBOW_SIGN_H
