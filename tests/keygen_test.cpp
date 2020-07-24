//
// Created by Florian Caullery on 7/24/20.
//

#include "gtest/gtest.h"
#include "../src/lib/error_codes.h"


extern "C" {
#include "../src/main/timing-functions.h"
#include "../src/lib/keygen.h"
#include "../src/lib/utils_prng.h"
#include "../src/lib/utils_hash.h"
}
extern "C" {
#include "../src/lib/gf16.h"
}

TEST(keygen, generation_of_s) {
    matrix_s_t s;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_matrix_s(s, &prng);

    uint64_t test_value = 0x01u;
    EXPECT_EQ(return_value, SUCCESS);
    for (int i = 0; i < O1; i++) {
        EXPECT_EQ(s[i].c, test_value);
        EXPECT_EQ(s[i].x, 0);
        EXPECT_EQ(s[i].y, 0);
        EXPECT_EQ(s[i].y_x, 0);
        test_value <<= 1u;
    }
    uint64_t mask_for_upper_value = 0xFFFFFFFF00000000u;
    for (int i = O1; i < O1 + O2; i++) {
        EXPECT_EQ(s[i].c & mask_for_upper_value, test_value);
        EXPECT_EQ(s[i].x & mask_for_upper_value, 0);
        EXPECT_EQ(s[i].y & mask_for_upper_value, 0);
        EXPECT_EQ(s[i].y_x & mask_for_upper_value, 0);
        test_value <<= 1u;
    }
}

TEST(keygen, generation_of_t) {
    matrix_t_t t;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_matrix_t(t, &prng);

    uint64_t test_value = 0x01u;
    EXPECT_EQ(return_value, SUCCESS);
    for (int i = 0; i < V1; i++) {
        EXPECT_EQ(t[0][i].c, test_value);
        EXPECT_EQ(t[0][i].x, 0);
        EXPECT_EQ(t[0][i].y, 0);
        EXPECT_EQ(t[0][i].y_x, 0);
        EXPECT_EQ(t[1][i].c, 0);
        EXPECT_EQ(t[1][i].x, 0);
        EXPECT_EQ(t[1][i].y, 0);
        EXPECT_EQ(t[1][i].y_x, 0);
        test_value <<= 1u;
    }

    uint64_t mask_for_upper_value = 0xFFFFFFFF00000000u;
    for (int i = V1; i < V1 + O1; i++) {
        EXPECT_EQ(t[0][i].c & mask_for_upper_value, test_value);
        EXPECT_EQ(t[0][i].x & mask_for_upper_value, 0);
        EXPECT_EQ(t[0][i].y & mask_for_upper_value, 0);
        EXPECT_EQ(t[0][i].y_x & mask_for_upper_value, 0);
        EXPECT_EQ(t[1][i].c, 0);
        EXPECT_EQ(t[1][i].x, 0);
        EXPECT_EQ(t[1][i].y, 0);
        EXPECT_EQ(t[1][i].y_x, 0);
        test_value <<= 1u;
    }
    test_value = 1u;
    for (int i = V1 + O1; i < V1 + O1 + O2; i++) {
        EXPECT_EQ(t[1][i].c, test_value);
        EXPECT_EQ(t[1][i].x, 0);
        EXPECT_EQ(t[1][i].y, 0);
        EXPECT_EQ(t[1][i].y_x, 0);
        test_value <<= 1u;
    }
}

TEST(keygen, generate_random_matrices_f) {
    matrix_fi_t f[O1 + O2];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_matrices_f(f, &prng);


    uint64_t test_value = 0xFFFFFFFFlu;
    uint64_t i, j;
    uint64_t accumulator_1 = 0;
    uint64_t accumulator_2 = 0;
    for (i = 0; i < O1; i++) {
        for (j = 0; j < O1 - 1; j++) {
            EXPECT_EQ(f[i][0][j].c & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][0][j].y & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][0][j].x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][0][j].y_x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][1][j].c, 0);
            EXPECT_EQ(f[i][1][j].y, 0);
            EXPECT_EQ(f[i][1][j].x, 0);
            EXPECT_EQ(f[i][1][j].y_x, 0);
            accumulator_1 |= f[i][0][j].c | f[i][0][j].y | f[i][0][j].x | f[i][0][j].y_x;
            accumulator_2 |= f[i][1][j].c | f[i][1][j].y | f[i][1][j].x | f[i][1][j].y_x;
        }
        for (j = O1; j < O1 + O2; j++) {
            EXPECT_EQ(f[i][0][j].c & ~test_value, 0);
            EXPECT_EQ(f[i][0][j].y & ~test_value, 0);
            EXPECT_EQ(f[i][0][j].x & ~test_value, 0);
            EXPECT_EQ(f[i][0][j].y_x & ~test_value, 0);
            EXPECT_EQ(f[i][1][j].c, 0);
            EXPECT_EQ(f[i][1][j].y, 0);
            EXPECT_EQ(f[i][1][j].x, 0);
            EXPECT_EQ(f[i][1][j].y_x, 0);
            accumulator_1 |= f[i][0][j].c | f[i][0][j].y | f[i][0][j].x | f[i][0][j].y_x;
            accumulator_2 |= f[i][1][j].c | f[i][1][j].y | f[i][1][j].x | f[i][1][j].y_x;
        }
        for (j = O1 + O2; j < O1 + O2 + V1; j++) {
            EXPECT_EQ(f[i][0][j].c, 0);
            EXPECT_EQ(f[i][0][j].y, 0);
            EXPECT_EQ(f[i][0][j].x, 0);
            EXPECT_EQ(f[i][0][j].y_x, 0);
            EXPECT_EQ(f[i][1][j].c, 0);
            EXPECT_EQ(f[i][1][j].y, 0);
            EXPECT_EQ(f[i][1][j].x, 0);
            EXPECT_EQ(f[i][1][j].y_x, 0);
            accumulator_1 |= f[i][0][j].c | f[i][0][j].y | f[i][0][j].x | f[i][0][j].y_x;
            accumulator_2 |= f[i][1][j].c | f[i][1][j].y | f[i][1][j].x | f[i][1][j].y_x;
        }
    }
    EXPECT_NE(accumulator_1, 0);
    EXPECT_EQ(accumulator_2, 0);
    accumulator_1 = 0;
    accumulator_2 = 0;

    for (i = V1; i < O1 + O2; i++) {
        for (j = 0; j < O1 + O2 - 1; j++) {
            EXPECT_EQ(f[i][0][j].c & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][0][j].y & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][0][j].x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][0][j].y_x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(f[i][1][j].c, 0);
            EXPECT_EQ(f[i][1][j].y, 0);
            EXPECT_EQ(f[i][1][j].x, 0);
            EXPECT_EQ(f[i][1][j].y_x, 0);
            accumulator_1 |= f[i][0][j].c | f[i][0][j].y | f[i][0][j].x | f[i][0][j].y_x;
            accumulator_2 |= f[i][1][j].c | f[i][1][j].y | f[i][1][j].x | f[i][1][j].y_x;
        }
        for (j = O1 + O2; j < O1 + O2 + V1; j++) {
            EXPECT_EQ(f[i][1][j].c, 0);
            EXPECT_EQ(f[i][1][j].y, 0);
            EXPECT_EQ(f[i][1][j].x, 0);
            EXPECT_EQ(f[i][1][j].y_x, 0);
            accumulator_1 |= f[i][0][j].c | f[i][0][j].y | f[i][0][j].x | f[i][0][j].y_x;
            accumulator_2 |= f[i][1][j].c | f[i][1][j].y | f[i][1][j].x | f[i][1][j].y_x;
        }
    }
    EXPECT_NE(accumulator_1, 0);
    EXPECT_EQ(accumulator_2, 0);
}

TEST(keygen, transpose_and_add_32x32_gf16_matrices) {

    uint64_t value = 1lu;
    unsigned int i;
    matrix_s_t s1, s1_plus_s1t;
    memset(s1, 0, sizeof(matrix_s_t));
    memset(s1_plus_s1t, 0, sizeof(matrix_s_t));
    for (i = 0; i < 32; i++) {
        s1[i].c = value;
        s1[i].y = value;
        s1[i].x = value;
        s1[i].y_x = value;
        value <<= 1u;
        value |= 1u;
    }
    transpose_and_add_32x32_gf16_matrices(s1_plus_s1t, s1, s1);
    for (i = 0; i < 32; i++) {
        EXPECT_EQ(s1_plus_s1t[i].c, 0xFFFFFFFF ^ (1u << i));
        EXPECT_EQ(s1_plus_s1t[i].y, 0xFFFFFFFF ^ (1u << i));
        EXPECT_EQ(s1_plus_s1t[i].x, 0xFFFFFFFF ^ (1u << i));
        EXPECT_EQ(s1_plus_s1t[i].y_x, 0xFFFFFFFF ^ (1u << i));
    }
}