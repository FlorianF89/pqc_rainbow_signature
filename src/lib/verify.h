//
// Created by Florian Caullery on 8/10/20.
//

#ifndef PQC_RAINBOW_VERIFY_H
#define PQC_RAINBOW_VERIFY_H


#include "gf16.h"
#include "keygen.h"

void evaluate_quadratic_polynomials(bitsliced_gf16_t *evaluations,
                                    bitsliced_quadratic_polynomials_t *f,
                                    bitsliced_gf16_t *x0_x63, bitsliced_gf16_t *x64_x95);
#endif //PQC_RAINBOW_VERIFY_H
