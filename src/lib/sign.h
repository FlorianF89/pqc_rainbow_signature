//
// Created by Florian Caullery on 8/2/20.
//

#ifndef PQC_RAINBOW_SIGN_H
#define PQC_RAINBOW_SIGN_H

#include "gf16.h"
#include "utils_prng.h"
#include "keygen.h"

void evaluate_quadratic_polynomials_at_x0_x31(bitsliced_gf16_t *evaluations, bitsliced_quadratic_polynomials_t *f,
                                              bitsliced_gf16_t *x0_x31);

void evaluate_quadratic_polynomials_at_x0_x63(bitsliced_gf16_t *evaluations, bitsliced_quadratic_polynomials_t *f,
                                              bitsliced_gf16_t *x0_x63);


void evaluate_32_quadratic_polynomials_at_x0_x63(bitsliced_gf16_t *evaluation, bitsliced_gf16_t f[32][64],
                                                 bitsliced_gf16_t *x0_x63);

int solve_32x32_gf16_system(bitsliced_gf16_t *solution, bitsliced_gf16_t equations_coefficients[32],
                            bitsliced_gf16_t *linear_coefficients);

void evaluate_quadratic_polynomials_of_second_layer_at_x0_x95(bitsliced_gf16_t *evaluations,
                                                              bitsliced_quadratic_polynomials_t *f,
                                                              bitsliced_gf16_t *x0_x63, bitsliced_gf16_t *x64_x95);

void find_preimage_of_x0_x31_by_32_polynomials_of_first_layer(bitsliced_gf16_t *preimages,
                                                              bitsliced_gf16_t *evaluation_in_x0_x31,
                                                              bitsliced_quadratic_polynomials_t *f,
                                                              bitsliced_gf16_t *x0_x31, prng_t *prng);

int find_preimage_of_x64_x96_by_32_polynomials_of_second_layer(bitsliced_gf16_t *y64_y96,
                                                               bitsliced_quadratic_polynomials_t *f,
                                                               bitsliced_gf16_t *y0_y64, bitsliced_gf16_t *x0_x64,
                                                               bitsliced_gf16_t *evaluation_in_x0_x31);


#endif //PQC_RAINBOW_SIGN_H
