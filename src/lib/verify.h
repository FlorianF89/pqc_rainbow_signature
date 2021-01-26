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

int rainbow_verify(bitsliced_gf16_t *signature, public_key_t *public_key, uint8_t *message, uint64_t message_length);

#endif //PQC_RAINBOW_VERIFY_H
