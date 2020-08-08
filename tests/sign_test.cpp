//
// Created by Florian Caullery on 8/2/20.
//


#include "gtest/gtest.h"
#include "../src/lib/error_codes.h"
#include <x86intrin.h>


extern "C" {
#include "../src/main/timing-functions.h"
#include "../src/lib/keygen.h"
#include "../src/lib/utils_prng.h"
#include "../src/lib/utils_hash.h"
#include "../src/lib/rng.h"
#include "../src/lib/sign.h"
#include "../src/lib/gf16.h"
}

TEST(sign_tests, evaluate_quadratic_polynomial_at_x0_x31) {
    bitsliced_gf16_t x0_x31, evaluation;
    memset(&evaluation, 0xFF, sizeof(bitsliced_gf16_t));
    x0_x31.c = 0xFFFFFFFFlu;
    x0_x31.y = 0;
    x0_x31.x = 0;
    x0_x31.y_x = 0;

    bitsliced_gf16_t f[32];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_32x32_gf16_upper_triangular_matrix(f, &prng);
    ASSERT_EQ(return_value, SUCCESS);
    int i, j;
    bitsliced_gf16_t accumulator = f[0];
    for (i = 1; i < 32; i++) {
        bitsliced_addition(&accumulator, &accumulator, &f[i]);
    }
    bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&accumulator, &accumulator, 0);
    evaluate_quadratic_polynomial_at_x0_x31(&evaluation, 0, f, &x0_x31);
    bitsliced_addition(&evaluation, &evaluation, &accumulator);
    EXPECT_EQ(gf16_is_zero(evaluation, 0), 1);
}

TEST(sign_tests, evaluate_32_quadratic_polynomials_at_x0_x31) {
    bitsliced_gf16_t x0_x31, evaluation;
    memset(&evaluation, 0xFF, sizeof(bitsliced_gf16_t));
    x0_x31.c = 0xFFFFFFFFlu;
    x0_x31.y = 0;
    x0_x31.x = 0;
    x0_x31.y_x = 0;

    bitsliced_gf16_t f[32][32];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    int return_value = 0;
    for (i = 0; i < 32; i++) {
        return_value += generate_random_32x32_gf16_upper_triangular_matrix(f[i], &prng);
    }
    ASSERT_EQ(return_value, SUCCESS);
    bitsliced_gf16_t accumulator;
    memset(&accumulator, 0x00, sizeof(bitsliced_gf16_t));
    for (i = 0; i < 32; i++) {
        bitsliced_gf16_t tmp = f[i][0];
        for (j = 1; j < 32; j++) {
            bitsliced_addition(&tmp, &tmp, &f[i][j]);
        }
        bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, i);
        bitsliced_addition(&accumulator, &tmp, &accumulator);
    }
    evaluate_32_quadratic_polynomials_at_x0_x31(&evaluation, f, &x0_x31);
    bitsliced_addition(&evaluation, &evaluation, &accumulator);
    EXPECT_EQ(gf16_is_zero(evaluation, 0), 1);
}

TEST(sign_tests, evaluate_32_quadratic_polynomials_at_x0_x63) {
    bitsliced_gf16_t x0_x31, evaluation;
    memset(&evaluation, 0xFF, sizeof(bitsliced_gf16_t));
    x0_x31.c = 0xFFFFFFFFFFFFFFFFlu;
    x0_x31.y = 0;
    x0_x31.x = 0;
    x0_x31.y_x = 0;

    bitsliced_gf16_t f[32][64];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    int return_value = 0;
    for (i = 0; i < 32; i++) {
        return_value += generate_random_32x32_gf16_upper_triangular_matrix(f[i], &prng);
        return_value += generate_random_32x32_gf16_matrix(&f[i][32], &prng);
    }
    ASSERT_EQ(return_value, SUCCESS);
    bitsliced_gf16_t accumulator;
    memset(&accumulator, 0x00, sizeof(bitsliced_gf16_t));
    for (i = 0; i < 32; i++) {
        bitsliced_gf16_t tmp = f[i][0];
        for (j = 1; j < 64; j++) {
            bitsliced_addition(&tmp, &tmp, &f[i][j]);
        }
        bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, i);
        bitsliced_addition(&accumulator, &tmp, &accumulator);
    }
    evaluate_32_quadratic_polynomials_at_x0_x63(&evaluation, f, &x0_x31);
    bitsliced_addition(&evaluation, &evaluation, &accumulator);
    EXPECT_EQ(gf16_is_zero(evaluation, 0), 1);
}

TEST(sign_tests, solve_32x32_gf16_system) {


    int i;
    int has_solution = 0;
    bitsliced_gf16_t linear_coefficients, x0_x31_prime, solution;
    bitsliced_gf16_t f[32];
    bitsliced_gf16_t f_copy[32];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];

    memset(&x0_x31_prime, 0xff, sizeof(bitsliced_gf16_t));
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    clock_t t, total = 0;
    int iterations = 0;
    while (has_solution == 0) {

        prng_gen(&prng, (uint8_t *) &linear_coefficients, sizeof(bitsliced_gf16_t));
        linear_coefficients.c &= 0xFFFFFFFFlu;
        linear_coefficients.y &= 0xFFFFFFFFlu;
        linear_coefficients.x &= 0xFFFFFFFFlu;
        linear_coefficients.y_x &= 0xFFFFFFFFlu;

        int return_value = 0;
        for (i = 0; i < 32; i++) {
            return_value += generate_random_32x32_gf16_matrix(f, &prng);
        }
        memcpy(f_copy, f, sizeof(f));
        t = clock();
        has_solution = solve_32x32_gf16_system(&solution, f, &linear_coefficients);
        total += clock() - t;
        iterations++;
    }
    printf("solve: %lucc\n", total / iterations);
    bitsliced_gf16_t g[32], verif[32];
    for (i = 0; i < 32; i++) {
        copy_gf16(&g[i], &solution);
    }
    multiply_32x32_gf16_matrices(verif, f_copy, g);
    bitsliced_addition(&g[0], &verif[0], &linear_coefficients);
    for (i = 0; i < 32; i++) {
        EXPECT_EQ(gf16_is_zero(g[0], i), 1);
    }
}

TEST(sign_tests, find_preimage_of_x0_x31_by_32_polynomials_in_64_variables) {
    bitsliced_gf16_t x0_x31, x0_x31_prime, y0_y63;
    bitsliced_gf16_t f[32][64];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];

    memset(&x0_x31_prime, 0xff, sizeof(bitsliced_gf16_t));
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t *) &x0_x31, sizeof(bitsliced_gf16_t));
    x0_x31.c &= 0xFFFFFFFFlu;
    x0_x31.y &= 0xFFFFFFFFlu;
    x0_x31.x &= 0xFFFFFFFFlu;
    x0_x31.y_x &= 0xFFFFFFFFlu;

    int i, j;
    int return_value = 0;
    for (i = 0; i < 32; i++) {
        return_value += generate_random_32x32_gf16_upper_triangular_matrix(f[i], &prng);
        return_value += generate_random_32x32_gf16_matrix(&f[i][32], &prng);
    }

//    clock_t start, end, total = 0;
    unsigned int dummy;
    unsigned long long t1 = __rdtscp(&dummy);
    for (i = 0; i < 2; i++) {
        find_preimage_of_x0_x31_by_32_polynomials_in_64_variables(&y0_y63, f, &x0_x31, &prng);
        find_preimage_of_x0_x31_by_32_polynomials_in_64_variables(&y0_y63, f, &x0_x31, &prng);
    }
    unsigned long long t2 = __rdtscp(&dummy);
    std::cout << "Time: " << t2 - t1 << std::endl;
    evaluate_32_quadratic_polynomials_at_x0_x63(&x0_x31_prime, f, &y0_y63);
    bitsliced_addition(&x0_x31_prime, &x0_x31_prime, &x0_x31);
    EXPECT_EQ(gf16_is_zero(x0_x31_prime, 0), 1);
}

