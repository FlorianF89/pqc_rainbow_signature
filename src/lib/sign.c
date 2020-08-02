//
// Created by Florian Caullery on 8/2/20.
//

#include <memory.h>
#include "sign.h"
#include "keygen.h"


void evaluate_quadratic_polynomial_at_x0_x31(bitsliced_gf16_t *evaluation, uint8_t position_to_put_evaluation,
                                             bitsliced_gf16_t *f, bitsliced_gf16_t *x0_x31) {
    bitsliced_gf16_t f_transposed[32];
    bitsliced_gf16_t tmp;
    transpose_32x32_gf16_matrix(f_transposed, f);
    memset(evaluation, 0, sizeof(bitsliced_gf16_t));
    uint32_t i;
    for (i = 0; i < 32; i++) {
        bitsliced_multiplication(&tmp, &f_transposed[i], x0_x31);
        bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, i);
        bitsliced_addition(evaluation, evaluation, &tmp);
    }
    bitsliced_multiplication(&tmp, evaluation, x0_x31);
    bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(evaluation, &tmp,
                                                                            position_to_put_evaluation);
}