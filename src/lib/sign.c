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


void evaluate_32_quadratic_polynomials_at_x0_x31(bitsliced_gf16_t *evaluation, bitsliced_gf16_t f[32][32],
                                                 bitsliced_gf16_t *x0_x31) {

    uint32_t i;
    memset(evaluation, 0, sizeof(bitsliced_gf16_t));
    for (i = 0; i < 32; i++) {
        bitsliced_gf16_t fi_transposed[32];
        bitsliced_gf16_t tmp, accumulator;
        transpose_32x32_gf16_matrix(fi_transposed, f[i]);
        memset(&accumulator, 0, sizeof(bitsliced_gf16_t));
        uint32_t j;
        for (j = 0; j < 32; j++) {
            bitsliced_multiplication(&tmp, &fi_transposed[j], x0_x31);
            bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, j);
            bitsliced_addition(&accumulator, &accumulator, &tmp);
        }
        bitsliced_multiplication(&tmp, &accumulator, x0_x31);
        bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&accumulator, &tmp, i);
        bitsliced_addition(evaluation, &accumulator, evaluation);
    }
}


void evaluate_32_quadratic_polynomials_at_x0_x63(bitsliced_gf16_t *evaluation, bitsliced_gf16_t f[32][64],
                                                 bitsliced_gf16_t *x0_x31) {

    uint32_t i;
    for (i = 0; i < 32; i++) {
        bitsliced_gf16_t fi_transposed[32];
        bitsliced_gf16_t tmp;
        transpose_32x32_gf16_matrix(fi_transposed, f[i]);
        memset(evaluation, 0, sizeof(bitsliced_gf16_t));
        uint32_t j;
        for (j = 0; j < 32; j++) {
            bitsliced_multiplication(&tmp, &fi_transposed[j], x0_x31);
            bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, j);
            bitsliced_addition(evaluation, evaluation, &tmp);
        }
        bitsliced_multiplication(&tmp, evaluation, x0_x31);
        bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(evaluation, &tmp, i);
    }
}

void find_preimage_of_x0_x31_by_32_polynomials_in_64_variables(bitsliced_gf16_t *preimage, bitsliced_gf16_t f[32][32],
                                                               bitsliced_gf16_t *x0_x31) {

}

