//
// Created by Florian Caullery on 8/2/20.
//

#ifndef PQC_RAINBOW_SIGN_H
#define PQC_RAINBOW_SIGN_H

#include "gf16.h"

void evaluate_quadratic_polynomial_at_x0_x31(bitsliced_gf16_t *evaluation, uint8_t position_to_put_evaluation,
                                             bitsliced_gf16_t *f, bitsliced_gf16_t *x0_x31);

#endif //PQC_RAINBOW_SIGN_H
