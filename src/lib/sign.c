//
// Created by Florian Caullery on 8/2/20.
//

#include <memory.h>
#include "sign.h"
#include "keygen.h"

void evaluate_quadratic_polynomials_at_x0_x31(bitsliced_gf16_t *evaluations, bitsliced_quadratic_polynomials_t *f,
                                              bitsliced_gf16_t *x0_x31) {
    unsigned int i, j;
    bitsliced_gf16_t x0_to_x31_expanded[32];
    for (i = 0; i < 32; i++) {
        x0_to_x31_expanded[i].c = -1 * ((x0_x31->c >> i) & 0x01u);
        x0_to_x31_expanded[i].y = -1 * ((x0_x31->y >> i) & 0x01u);
        x0_to_x31_expanded[i].x = -1 * ((x0_x31->x >> i) & 0x01u);
        x0_to_x31_expanded[i].y_x = -1 * ((x0_x31->y_x >> i) & 0x01u);
    }

    evaluations->c = 0;
    evaluations->y = 0;
    evaluations->x = 0;
    evaluations->y_x = 0;

    bitsliced_gf16_t tmp, tmp1;
    for (i = 0; i < 32; i++) {
        for (j = i; j < 32; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            bitsliced_multiplication(&tmp, &x0_to_x31_expanded[i], &x0_to_x31_expanded[j]);
            bitsliced_multiplication(&tmp1, &tmp, &f->coefficients[position]);
            bitsliced_addition(evaluations, evaluations, &tmp1);
        }
    }
}

void evaluate_quadratic_polynomials_at_x0_x63(bitsliced_gf16_t *evaluations, bitsliced_quadratic_polynomials_t *f,
                                              bitsliced_gf16_t *x0_x64) {
    unsigned int i, j;
    bitsliced_gf16_t x0_to_x63_expanded[64];
    for (i = 0; i < 64; i++) {
        x0_to_x63_expanded[i].c = -1 * ((x0_x64->c >> i) & 0x01u);
        x0_to_x63_expanded[i].y = -1 * ((x0_x64->y >> i) & 0x01u);
        x0_to_x63_expanded[i].x = -1 * ((x0_x64->x >> i) & 0x01u);
        x0_to_x63_expanded[i].y_x = -1 * ((x0_x64->y_x >> i) & 0x01u);
    }

    evaluations->c = 0;
    evaluations->y = 0;
    evaluations->x = 0;
    evaluations->y_x = 0;

    bitsliced_gf16_t tmp, tmp1;
    for (i = 0; i < 64; i++) {
        for (j = i; j < 64; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            bitsliced_multiplication(&tmp, &x0_to_x63_expanded[i], &x0_to_x63_expanded[j]);
            bitsliced_multiplication(&tmp1, &tmp, &f->coefficients[position]);
            bitsliced_addition(evaluations, evaluations, &tmp1);
        }
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
        printf("l%02d: ", i);
        for (j = 0; j < 32; j++) {
            print_bitsliced(m[j], i);
        }
        print_bitsliced(solution, i);
        putchar(10);
    }
    putchar(10);
}

static bitsliced_gf16_t swap_elements_i_and_k(bitsliced_gf16_t in, unsigned int i, unsigned int j) {
    bitsliced_gf16_t r;
    uint64_t b = in.c;
    uint64_t x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r.c = b ^ ((x << i) | (x << j));

    b = in.y;
    x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r.y = b ^ ((x << i) | (x << j));

    b = in.x;
    x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r.x = b ^ ((x << i) | (x << j));

    b = in.y_x;
    x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r.y_x = b ^ ((x << i) | (x << j));

    return r;
}

//todo Does not work if system is undetermined
int solve_32x32_gf16_system(bitsliced_gf16_t *solution, bitsliced_gf16_t equations_coefficients[32],
                            bitsliced_gf16_t *linear_coefficients) {
    copy_gf16(solution, linear_coefficients);
    unsigned int i, j, k;
    for (i = 0; i < 31; i++) {
        for (j = i + 1; j < 31; j++) {
            if (gf16_is_zero(equations_coefficients[i], i)) {
//                printf("\n swapping lines %d and %d. Before:\n", i, j);
//                print_32x32_gf16_system(equations_coefficients, *solution);
                for (k = i; k < 32; k++)
                    equations_coefficients[k] = swap_elements_i_and_k(equations_coefficients[k], i, j);
                *solution = swap_elements_i_and_k(*solution, i, j);
//                printf("\n swapping lines %d and %d. After:\n", i, j);
//                print_32x32_gf16_system(equations_coefficients, *solution);
            } else {
                break;
            }
        }
        if (j == 32) {
            return 0;
        }
        bitsliced_gf16_t tmp = extract_then_expand(equations_coefficients[i], i, i);
        bitsliced_gf16_t tmp1;
        bitsliced_inversion(&tmp1, &tmp);
        bitsliced_multiplication(&tmp, &tmp1, &equations_coefficients[i]);
        for (j = i + 1u; j < 32; ++j) {
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
    for (i = 31; i > 0; i--) {
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

void find_preimage_of_x0_x31_by_32_polynomials_of_first_layer(bitsliced_gf16_t *preimages,
                                                              bitsliced_quadratic_polynomials_t *f,
                                                              bitsliced_gf16_t *x0_x31, prng_t *prng) {

    uint8_t i, j;
    bitsliced_gf16_t tmp, tmp2, accumulator;
    bitsliced_gf16_t linear_system[32];
    memset(linear_system, 0x00, sizeof(linear_system));
    int found = 0;
    int attempt = 0;
    while (found == 0 && attempt < 100000) {
        prng_gen(prng, (uint8_t *) preimages, sizeof(bitsliced_gf16_t));
        evaluate_quadratic_polynomials_at_x0_x31(&accumulator, f, preimages);
        bitsliced_addition(&accumulator, x0_x31, &accumulator);
        bitsliced_gf16_t preimages_0_31_expanded[32];
        for (i = 0; i < 32; i++) {
            preimages_0_31_expanded[i].c = -1 * ((preimages->c >> i) & 0x01u);
            preimages_0_31_expanded[i].y = -1 * ((preimages->y >> i) & 0x01u);
            preimages_0_31_expanded[i].x = -1 * ((preimages->x >> i) & 0x01u);
            preimages_0_31_expanded[i].y_x = -1 * ((preimages->y_x >> i) & 0x01u);
        }
        for (i = 0; i < 32; i++) {
            for (j = 32; j < 64; j++) {
                int position = i * N - ((i + 1) * i / 2) + j;
                bitsliced_multiplication(&tmp, &preimages_0_31_expanded[i], &f->coefficients[position]);
                bitsliced_addition(&linear_system[j - 32], &linear_system[j - 32], &tmp);
            }
        }
        found = solve_32x32_gf16_system(&tmp2, linear_system, &accumulator);
        if (found) {
            preimages->c = (preimages->c & 0xFFFFFFFF) | (tmp2.c << 32u);
            preimages->y = (preimages->y & 0xFFFFFFFF) | (tmp2.y << 32u);
            preimages->x = (preimages->x & 0xFFFFFFFF) | (tmp2.x << 32u);
            preimages->y_x = (preimages->y_x & 0xFFFFFFFF) | (tmp2.y_x << 32u);
        }
        attempt++;
    }
}

