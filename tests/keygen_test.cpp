//
// Created by Florian Caullery on 7/24/20.
//

#include <x86intrin.h>
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
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x |
                    private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x |
                    private_key.f[i][1][j].y_x;
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
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x |
                    private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x |
                    private_key.f[i][1][j].y_x;
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
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x |
                    private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x |
                    private_key.f[i][1][j].y_x;
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
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x |
                    private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x |
                    private_key.f[i][1][j].y_x;
        }
        for (j = O1 + O2; j < O1 + O2 + V1; j++) {
            EXPECT_EQ(private_key.f[i][1][j].c, 0);
            EXPECT_EQ(private_key.f[i][1][j].y, 0);
            EXPECT_EQ(private_key.f[i][1][j].x, 0);
            EXPECT_EQ(private_key.f[i][1][j].y_x, 0);
            accumulator_1 |=
                    private_key.f[i][0][j].c | private_key.f[i][0][j].y | private_key.f[i][0][j].x |
                    private_key.f[i][0][j].y_x;
            accumulator_2 |=
                    private_key.f[i][1][j].c | private_key.f[i][1][j].y | private_key.f[i][1][j].x |
                    private_key.f[i][1][j].y_x;
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
    int total_iterations = 100;
    uint32_t dummy;
    unsigned long long t1, t2, total = 0;
    for (j = 0; j < total_iterations; j++) {
        for (i = 0; i < 32; i++) {
            prng_gen(&prng, (unsigned char *) &s[i].c, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &s[i].y, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &s[i].x, sizeof(uint32_t));
            prng_gen(&prng, (unsigned char *) &s[i].y_x, sizeof(uint32_t));
        }
        t1 = __rdtscp(&dummy);
        gaussian_elimination_for_32x32_gf16_matrix(s);
        t2 = __rdtscp(&dummy);
        total += t2 - t1;
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
    int i, j, k, l;
    clock_t t;
    clock_t total = 0;
    int total_iterations = 1;
    for (j = 0; j < total_iterations; j++) {
        t = clock();
        generate_private_key(&private_key, &prng);
        total += clock() - t;
        uint64_t s_accumulator = 0;
        uint64_t first_layer_accumulator = 0;
        uint64_t second_layer_accumulator = 0;
        //testing first layer
        for (i = 0; i < N; i++) {
            for (k = i; k < N; k++) {
                int position = i * N - ((i + 1) * i / 2) + k;
                if (i >= O1 || k >= O1 + O2) {
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].c & 0xFFFFFFFF, 0);
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].y & 0xFFFFFFFF, 0);
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].x & 0xFFFFFFFF, 0);
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].y_x & 0xFFFFFFFF, 0);
                } else {
                    first_layer_accumulator |= private_key.mq_polynomials.coefficients[position].c & 0xFFFFFFFF;
                    first_layer_accumulator |= private_key.mq_polynomials.coefficients[position].y & 0xFFFFFFFF;
                    first_layer_accumulator |= private_key.mq_polynomials.coefficients[position].x & 0xFFFFFFFF;
                    first_layer_accumulator |= private_key.mq_polynomials.coefficients[position].y_x & 0xFFFFFFFF;
                }
            }
        }
        //testing second layer
        for (i = 0; i < N; i++) {
            for (k = i; k < N; k++) {
                int position = i * N - ((i + 1) * i / 2) + k;
                if (i >= O1 + O2 || k >= O1 + O2) {
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].c & 0xFFFFFFFF00000000, 0);
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].y & 0xFFFFFFFF00000000, 0);
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].x & 0xFFFFFFFF00000000, 0);
                    EXPECT_EQ(private_key.mq_polynomials.coefficients[position].y_x & 0xFFFFFFFF00000000, 0);
                } else {
                    second_layer_accumulator |=
                            private_key.mq_polynomials.coefficients[position].c & 0xFFFFFFFF00000000;
                    second_layer_accumulator |=
                            private_key.mq_polynomials.coefficients[position].y & 0xFFFFFFFF00000000;
                    second_layer_accumulator |=
                            private_key.mq_polynomials.coefficients[position].x & 0xFFFFFFFF00000000;
                    second_layer_accumulator |=
                            private_key.mq_polynomials.coefficients[position].y_x & 0xFFFFFFFF00000000;
                }
            }
        }
        EXPECT_NE(first_layer_accumulator, 0);
        EXPECT_NE(second_layer_accumulator, 0);
    }
    printf("private key gen: %f \n", (double) total / total_iterations);
}

TEST(keygen, derive_public_key_from_private_key) {

    private_key_t private_key;
    memset(&private_key, 0x00, sizeof(private_key));
    public_key_t public_key;
    memset(&public_key, 0x00, sizeof(public_key));
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    clock_t t;
//    clock_t total = 0;
    int total_iterations = 10;
    unsigned int dummy;
    unsigned long long t1, t2, total = 0;
    for (j = 0; j < total_iterations; j++) {
//        t = clock();
        generate_private_key(&private_key, &prng);
        t1 = __rdtscp(&dummy);
        derive_public_key_from_private_key(&public_key, &private_key);
        t2 = __rdtscp(&dummy);
        total += t2 - t1;
//        total += clock() - t;
        uint64_t mq_accumulator = 0;
        uint64_t mp_accumulator = 0;
        for (i = 0; i < (N * (N + 1) / 2); i++) {
            mq_accumulator |= public_key.mq[i].c;
            mp_accumulator |= public_key.mp[i].c;
            mq_accumulator |= public_key.mq[i].y;
            mp_accumulator |= public_key.mp[i].y;
            mq_accumulator |= public_key.mq[i].x;
            mp_accumulator |= public_key.mp[i].x;
            mq_accumulator |= public_key.mq[i].y_x;
            mp_accumulator |= public_key.mp[i].y_x;
        }
        EXPECT_NE(mq_accumulator, 0);
        EXPECT_NE(mp_accumulator, 0);
    }
    printf("public key gen: %f \n", (double) total / total_iterations);
}

TEST(keygen, replace_variable_by_linear_combination_in_quadratic_polynomial) {

    bitsliced_quadratic_polynomials_t f, f_prime;
    memset(&f, 0x00, sizeof(bitsliced_quadratic_polynomials_t));
    f.coefficients[0].c = -1;
    bitsliced_gf16_t linear_combination;
    linear_combination.c = -1;
    linear_combination.y = 0;
    linear_combination.x = 0;
    linear_combination.y_x = 0;

    unsigned int dummy;
    unsigned long long t1, t2, total = 0;
    t1 = __rdtscp(&dummy);
    replace_variable_by_linear_combination_in_quadratic_polynomial(&f_prime, &f, 0, &linear_combination);
    t2 = __rdtscp(&dummy);
    total += t2 - t1;
    printf("replace variables: %f \n", (double) total * 32);
    int i;
    int next_square_coeff = 32 * N - (32 * 33) / 2 + 32;
    int square_coeff_left = 64;
    EXPECT_EQ(f_prime.coefficients[0].c, -1);
    EXPECT_EQ(f_prime.coefficients[0].y, 0);
    EXPECT_EQ(f_prime.coefficients[0].x, 0);
    EXPECT_EQ(f_prime.coefficients[0].y_x, 0);
    for (i = 1; i < N * (N + 1) / 2; i++) {
        if (i == next_square_coeff) {
            EXPECT_EQ(f_prime.coefficients[i].c, -1);
            EXPECT_EQ(f_prime.coefficients[i].y, 0);
            EXPECT_EQ(f_prime.coefficients[i].x, 0);
            EXPECT_EQ(f_prime.coefficients[i].y_x, 0);
            next_square_coeff += square_coeff_left;
            square_coeff_left--;
        } else {
            EXPECT_EQ(f_prime.coefficients[i].c, 0);
            EXPECT_EQ(f_prime.coefficients[i].y, 0);
            EXPECT_EQ(f_prime.coefficients[i].x, 0);
            EXPECT_EQ(f_prime.coefficients[i].y_x, 0);
        }
    }
}