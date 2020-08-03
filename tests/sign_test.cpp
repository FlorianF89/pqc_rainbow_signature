//
// Created by Florian Caullery on 8/2/20.
//


#include "gtest/gtest.h"
#include "../src/lib/error_codes.h"


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

    find_preimage_of_x0_x31_by_32_polynomials_in_64_variables(&y0_y63, f, &x0_x31);
    evaluate_32_quadratic_polynomials_at_x0_x63(&x0_x31_prime, f, &y0_y63);
    bitsliced_addition(&x0_x31_prime, &x0_x31_prime, &x0_x31);
    EXPECT_EQ(gf16_is_zero(x0_x31_prime, 0), 1);
}

