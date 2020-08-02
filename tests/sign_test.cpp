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