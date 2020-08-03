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
                                                 bitsliced_gf16_t *x0_x63) {

    uint32_t i;
    memset(evaluation, 0, sizeof(bitsliced_gf16_t));
    for (i = 0; i < 32; i++) {
        bitsliced_gf16_t fi_left_transposed[32];
        bitsliced_gf16_t fi_right_transposed[32];
        bitsliced_gf16_t tmp, accumulator;
        transpose_32x32_gf16_matrix(fi_left_transposed, f[i]);
        transpose_32x32_gf16_matrix(fi_right_transposed, &f[i][32]);
        uint32_t j;
        for (j = 0; j < 32; j++) {
            shift_left_gf16(&fi_right_transposed[j], &fi_right_transposed[j], 32);
            bitsliced_addition(&fi_left_transposed[j], &fi_left_transposed[j], &fi_right_transposed[j]);
        }
        memset(&accumulator, 0, sizeof(bitsliced_gf16_t));
        for (j = 0; j < 32; j++) {
            bitsliced_multiplication(&tmp, &fi_left_transposed[j], x0_x63);
            bitsliced_gf16_sum_64_first_elements_and_place_result_in_given_position(&tmp, &tmp, j);
            bitsliced_addition(&accumulator, &accumulator, &tmp);
        }
        bitsliced_multiplication(&tmp, &accumulator, x0_x63);
        bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&accumulator, &tmp, i);
        bitsliced_addition(evaluation, &accumulator, evaluation);
    }
}


static void print_bitsliced(const bitsliced_gf16_t in, unsigned int bit_position) {
    if ((in.c >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in.y >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in.x >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in.y_x >> bit_position) & 0x01u) {
        printf("1 ");
    } else {
        printf("0 ");
    }
}

void print_32x32_gf16_system(bitsliced_gf16_t m[32], bitsliced_gf16_t solution) {
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            print_bitsliced(m[j], i);
        }
        print_bitsliced(solution, i);
        putchar(10);
    }
    putchar(10);
}

int solve_32x32_gf16_system(bitsliced_gf16_t *solution, bitsliced_gf16_t equations_coefficients[32],
                            bitsliced_gf16_t *linear_coefficients) {
    copy_gf16(solution, linear_coefficients);
    print_32x32_gf16_system(equations_coefficients, *linear_coefficients);
    for (unsigned int i = 0; i < 31; i++) {
        bitsliced_gf16_t tmp = extract_then_expand(equations_coefficients[i], i, i);
        bitsliced_gf16_t tmp1;
        bitsliced_inversion(&tmp1, &tmp);
        bitsliced_multiplication(&tmp, &tmp1, &equations_coefficients[i]);
        for (unsigned int j = i + 1u; j < 32; ++j) {
            bitsliced_gf16_t tmp2 = extract_then_expand(equations_coefficients[j], i, i + 1);
            bitsliced_multiplication(&tmp1, &tmp, &tmp2);
            bitsliced_addition(&equations_coefficients[j], &tmp1, &equations_coefficients[j]);
        }
        bitsliced_gf16_t tmp2 = extract_then_expand(*solution, i, i + 1);
        bitsliced_multiplication(&tmp1, &tmp, &tmp2);
        bitsliced_addition(solution, &tmp1, solution);

        equations_coefficients[i].c &= (1u << (i + 1u)) - 1u;
        equations_coefficients[i].y &= (1u << (i + 1u)) - 1u;
        equations_coefficients[i].x &= (1u << (i + 1u)) - 1u;
        equations_coefficients[i].y_x &= (1u << (i + 1u)) - 1u;
    }
    print_32x32_gf16_system(equations_coefficients, *solution);

    for (unsigned int i = 31; i > 0; i--) {
        bitsliced_gf16_t tmp = extract_then_expand(equations_coefficients[i], i, 0);
        shift_right_gf16(&tmp, &tmp, 64 - i);
        bitsliced_gf16_t tmp1, tmp2;
        bitsliced_inversion(&tmp1, &tmp);
        bitsliced_multiplication(&tmp, &tmp1, &equations_coefficients[i]);
        tmp1 = extract_then_expand(*solution, i, 0);
        shift_right_gf16(&tmp1, &tmp1, 64 - i);
        bitsliced_multiplication(&tmp2, &tmp, &tmp1);
        bitsliced_addition(solution, &tmp2, solution);
    }
    bitsliced_gf16_t tmp, diagonal_elements;
    memset(&diagonal_elements, 0x00, sizeof(bitsliced_gf16_t));
    int has_solution = 0;
    for (unsigned int i = 0; i < 32; i++) {
        extract_one_gf16_element_and_place_it_in_given_position(&tmp, i, &equations_coefficients[i], i);
        has_solution |= gf16_is_zero(tmp, i);
        bitsliced_addition(&diagonal_elements, &diagonal_elements, &tmp);
    }
    bitsliced_inversion(&tmp, &diagonal_elements);
    copy_gf16(&diagonal_elements, solution);
    bitsliced_multiplication(solution, &tmp, &diagonal_elements);
    return has_solution == 0;
}

void find_preimage_of_x0_x31_by_32_polynomials_in_64_variables(bitsliced_gf16_t *preimage, bitsliced_gf16_t f[32][64],
                                                               bitsliced_gf16_t *x0_x31) {

}

