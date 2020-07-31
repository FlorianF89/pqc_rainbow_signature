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
#include "../src/lib/rng.h"
}
extern "C" {
#include "../src/lib/gf16.h"
}

TEST(keygen, generation_of_s) {
    private_key_t private_key;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_matrix_s(&private_key, &prng);

    uint64_t test_value = 0x01u;
    EXPECT_EQ(return_value, SUCCESS);
    for (int i = 0; i < O1; i++) {
        EXPECT_EQ(private_key.s[i].c, test_value);
        EXPECT_EQ(private_key.s[i].x, 0);
        EXPECT_EQ(private_key.s[i].y, 0);
        EXPECT_EQ(private_key.s[i].y_x, 0);
        test_value <<= 1u;
    }
    uint64_t mask_for_upper_value = 0xFFFFFFFF00000000u;
    for (int i = O1; i < O1 + O2; i++) {
        EXPECT_EQ(private_key.s[i].c & mask_for_upper_value, test_value);
        EXPECT_EQ(private_key.s[i].x & mask_for_upper_value, 0);
        EXPECT_EQ(private_key.s[i].y & mask_for_upper_value, 0);
        EXPECT_EQ(private_key.s[i].y_x & mask_for_upper_value, 0);
        test_value <<= 1u;
    }
}

TEST(keygen, generation_of_t) {
    private_key private_key;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_matrix_t(&private_key, &prng);

    uint64_t test_value = 0x01u;
    EXPECT_EQ(return_value, SUCCESS);
    for (int i = 0; i < V1; i++) {
        EXPECT_EQ(private_key.t[0][i].c, test_value);
        EXPECT_EQ(private_key.t[0][i].x, 0);
        EXPECT_EQ(private_key.t[0][i].y, 0);
        EXPECT_EQ(private_key.t[0][i].y_x, 0);
        EXPECT_EQ(private_key.t[1][i].c, 0);
        EXPECT_EQ(private_key.t[1][i].x, 0);
        EXPECT_EQ(private_key.t[1][i].y, 0);
        EXPECT_EQ(private_key.t[1][i].y_x, 0);
        test_value <<= 1u;
    }

    uint64_t mask_for_upper_value = 0xFFFFFFFF00000000u;
    for (int i = V1; i < V1 + O1; i++) {
        EXPECT_EQ(private_key.t[0][i].c & mask_for_upper_value, test_value);
        EXPECT_EQ(private_key.t[0][i].x & mask_for_upper_value, 0);
        EXPECT_EQ(private_key.t[0][i].y & mask_for_upper_value, 0);
        EXPECT_EQ(private_key.t[0][i].y_x & mask_for_upper_value, 0);
        EXPECT_EQ(private_key.t[1][i].c, 0);
        EXPECT_EQ(private_key.t[1][i].x, 0);
        EXPECT_EQ(private_key.t[1][i].y, 0);
        EXPECT_EQ(private_key.t[1][i].y_x, 0);
        test_value <<= 1u;
    }
    test_value = 1u;
    for (int i = V1 + O1; i < V1 + O1 + O2; i++) {
        EXPECT_EQ(private_key.t[1][i].c, test_value);
        EXPECT_EQ(private_key.t[1][i].x, 0);
        EXPECT_EQ(private_key.t[1][i].y, 0);
        EXPECT_EQ(private_key.t[1][i].y_x, 0);
        test_value <<= 1u;
    }
}

TEST(keygen, generate_random_matrices_f) {
    private_key_t private_key;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int return_value = generate_random_matrices_f(&private_key, &prng);
    EXPECT_EQ(return_value, SUCCESS);

    uint64_t test_value = 0xFFFFFFFFlu;
    uint64_t i, j;
    uint64_t accumulator_1 = 0;
    uint64_t accumulator_2 = 0;
    for (i = 0; i < O1; i++) {
        for (j = 0; j < O1 - 1; j++) {
            EXPECT_EQ(private_key.f[i][0][j].c & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][0][j].y & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][0][j].x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][0][j].y_x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][1][j].c, 0);
            EXPECT_EQ(private_key.f[i][1][j].y, 0);
            EXPECT_EQ(private_key.f[i][1][j].x, 0);
            EXPECT_EQ(private_key.f[i][1][j].y_x, 0);
            accumulator_1 |=
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x | private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x | private_key.f[i][1][j].y_x;
        }
        for (j = O1; j < O1 + O2; j++) {
            EXPECT_EQ(private_key.f[i][0][j].c & ~test_value, 0);
            EXPECT_EQ(private_key.f[i][0][j].y & ~test_value, 0);
            EXPECT_EQ(private_key.f[i][0][j].x & ~test_value, 0);
            EXPECT_EQ(private_key.f[i][0][j].y_x & ~test_value, 0);
            EXPECT_EQ(private_key.f[i][1][j].c, 0);
            EXPECT_EQ(private_key.f[i][1][j].y, 0);
            EXPECT_EQ(private_key.f[i][1][j].x, 0);
            EXPECT_EQ(private_key.f[i][1][j].y_x, 0);
            accumulator_1 |=
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x | private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x | private_key.f[i][1][j].y_x;
        }
        for (j = O1 + O2; j < O1 + O2 + V1; j++) {
            EXPECT_EQ(private_key.f[i][0][j].c, 0);
            EXPECT_EQ(private_key.f[i][0][j].y, 0);
            EXPECT_EQ(private_key.f[i][0][j].x, 0);
            EXPECT_EQ(private_key.f[i][0][j].y_x, 0);
            EXPECT_EQ(private_key.f[i][1][j].c, 0);
            EXPECT_EQ(private_key.f[i][1][j].y, 0);
            EXPECT_EQ(private_key.f[i][1][j].x, 0);
            EXPECT_EQ(private_key.f[i][1][j].y_x, 0);
            accumulator_1 |=
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x | private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x | private_key.f[i][1][j].y_x;
        }
    }
    EXPECT_NE(accumulator_1, 0);
    EXPECT_EQ(accumulator_2, 0);
    accumulator_1 = 0;
    accumulator_2 = 0;

    for (i = V1; i < O1 + O2; i++) {
        for (j = 0; j < O1 + O2 - 1; j++) {
            EXPECT_EQ(private_key.f[i][0][j].c & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][0][j].y & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][0][j].x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][0][j].y_x & ~((1lu << (j + 1)) - 1), 0);
            EXPECT_EQ(private_key.f[i][1][j].c, 0);
            EXPECT_EQ(private_key.f[i][1][j].y, 0);
            EXPECT_EQ(private_key.f[i][1][j].x, 0);
            EXPECT_EQ(private_key.f[i][1][j].y_x, 0);
            accumulator_1 |=
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x | private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x | private_key.f[i][1][j].y_x;
        }
        for (j = O1 + O2; j < O1 + O2 + V1; j++) {
            EXPECT_EQ(private_key.f[i][1][j].c, 0);
            EXPECT_EQ(private_key.f[i][1][j].y, 0);
            EXPECT_EQ(private_key.f[i][1][j].x, 0);
            EXPECT_EQ(private_key.f[i][1][j].y_x, 0);
            accumulator_1 |=
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x | private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x | private_key.f[i][1][j].y_x;
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

TEST(keygen, guassian_elimination_32x32_gf16_matrices) {

    bitsliced_gf16_t s[32];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    memset(s, 0x00, sizeof(s));
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    clock_t t;
    clock_t total = 0;
    int total_iterations = 10;
    for (j = 0; j < total_iterations; j++) {
        for (i = 0; i < 32; i++) {
            prng_gen(&prng, (unsigned char *) &s[i].c, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &s[i].y, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &s[i].x, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &s[i].y_x, sizeof(uint32_t));
        }
        t = clock();
        gaussian_elimination_for_32x32_gf16_matrix(s);
        total += clock() - t;
        uint64_t mask = 0xFFFFFFFE;
        for (i = 0; i < 32; i++) {
            EXPECT_EQ(s[i].c & mask, 0);
            EXPECT_EQ(s[i].y & mask, 0);
            EXPECT_EQ(s[i].x & mask, 0);
            EXPECT_EQ(s[i].y_x & mask, 0);
            mask <<= 1u;
        }
    }
    printf("%f \n", (double) total / total_iterations);
}


TEST(keygen, multiply_32x32_gf16_matrices) {

    bitsliced_gf16_t a[32], a_times_b[32], identity[32];
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    memset(a, 0x00, sizeof(a));
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    clock_t t;
    clock_t total = 0;
    int total_iterations = 10;
    set_32x32_gf16_matrix_to_identity(identity);
    for (j = 0; j < total_iterations; j++) {
        for (i = 0; i < 32; i++) {
            prng_gen(&prng, (unsigned char *) &a[i].c, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &a[i].y, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &a[i].x, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &a[i].y_x, sizeof(uint32_t));
        }
        t = clock();
        multiply_32x32_gf16_matrices(a_times_b, a, identity);
        total += clock() - t;
        for (i = 0; i < 32; i++) {
            EXPECT_EQ(a[i].c, a_times_b[i].c);
            EXPECT_EQ(a[i].y, a_times_b[i].y);
            EXPECT_EQ(a[i].x, a_times_b[i].x);
            EXPECT_EQ(a[i].y_x, a_times_b[i].y_x);
        }
    }
    printf("matrix mul: %f \n", (double) total / total_iterations);
}

TEST(keygen, generate_private_key) {

    private_key_t private_key;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    clock_t t;
    clock_t total = 0;
    int total_iterations = 10;
    for (j = 0; j < total_iterations; j++) {
        t = clock();
        generate_private_key(&private_key, &prng);
        total += clock() - t;
        uint64_t s_accumulator = 0;
        uint64_t t_accumulator = 0;
        uint64_t f_accumulator = 0;
        for (i = 0; i < 32; i++) {
            s_accumulator |= private_key.s[i + O1].c;
            s_accumulator |= private_key.s[i + O1].x;
            s_accumulator |= private_key.s[i + O1].y;
            s_accumulator |= private_key.s[i + O1].y_x;
            t_accumulator |= private_key.t[0][i + O1].c;
            t_accumulator |= private_key.t[0][i + O1].x;
            t_accumulator |= private_key.t[0][i + O1].y;
            t_accumulator |= private_key.t[0][i + O1].y_x;
            f_accumulator |= private_key.f[0][0][i + O1].c;
            f_accumulator |= private_key.f[0][0][i + O1].x;
            f_accumulator |= private_key.f[0][0][i + O1].y;
            f_accumulator |= private_key.f[0][0][i + O1].y_x;
        }
        EXPECT_NE(s_accumulator, 0);
        EXPECT_NE(t_accumulator, 0);
        EXPECT_NE(f_accumulator, 0);
    }
    printf("private key gen: %f \n", (double) total / total_iterations);
}

