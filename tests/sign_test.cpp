//
// Created by Florian Caullery on 8/2/20.
//


#include "gtest/gtest.h"
#include "../src/lib/error_codes.h"
#include <x86intrin.h>
#include <chrono>

template<class T>
__attribute__((always_inline)) inline void DoNotOptimize(const T &value) {
    asm volatile("" : "+m"(const_cast<T &>(value)));
}

extern "C" {
#include "../src/main/timing-functions.h"
#include "../src/lib/keygen.h"
#include "../src/lib/utils_prng.h"
#include "../src/lib/utils_hash.h"
#include "../src/lib/rng.h"
#include "../src/lib/sign.h"
#include "../src/lib/gf16.h"
}

auto time_sign() {

    private_key_t private_key;

    bitsliced_gf16_t signature[2];
    memset(signature, 0x00, sizeof(signature));

    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    uint8_t message[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    memset(message, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);

    int i;


    using Clock = std::chrono::high_resolution_clock;

    auto input = SECRET_KEY_SEED_BYTE_LENGTH;

    auto t1 = Clock::now();         // Statement 1
    DoNotOptimize(input);
    auto output = 0;
    for (i = 0; i < 20; i++) {
        output = rainbow_sign(signature, &private_key, message, &prng, input);
    }
    DoNotOptimize(output);
    auto t2 = Clock::now();         // Statement 3

    return t2 - t1;
}

TEST(sign_tests, evaluate_quadratic_polynomials_at_x0_x31) {
    bitsliced_gf16_t x0_x31, evaluation;
    memset(&evaluation, 0xFF, sizeof(bitsliced_gf16_t));
    x0_x31.c = 0xFFFFFFFFlu;
    x0_x31.y = 0;
    x0_x31.x = 0;
    x0_x31.y_x = 0;

    bitsliced_quadratic_polynomials_t f;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t *) f.coefficients, (N * (N + 1) / 2) * sizeof(bitsliced_gf16_t));
    int i, j;
    bitsliced_gf16_t accumulator;
    bitsliced_addition(&accumulator, &accumulator, &accumulator); //set acc to 0
    for (i = 0; i < 32; i++) {
        for (j = i; j < 32; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            bitsliced_addition(&accumulator, &accumulator, &f.coefficients[position]);
        }
    }
    evaluate_quadratic_polynomials_at_x0_x31(&evaluation, &f, &x0_x31);
    bitsliced_addition(&evaluation, &evaluation, &accumulator);
    EXPECT_EQ(evaluation.c, 0);
    EXPECT_EQ(evaluation.y, 0);
    EXPECT_EQ(evaluation.x, 0);
    EXPECT_EQ(evaluation.y_x, 0);
}

TEST(sign_tests, evaluate_quadratic_polynomials_at_x0_x63) {
    bitsliced_gf16_t x0_x63, evaluation, eval_x0_x31;
    memset(&evaluation, 0xFF, sizeof(bitsliced_gf16_t));
    x0_x63.c = -1;
    x0_x63.y = 0;
    x0_x63.x = 0;
    x0_x63.y_x = 0;

    bitsliced_quadratic_polynomials_t f;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t *) f.coefficients, (N * (N + 1) / 2) * sizeof(bitsliced_gf16_t));
    int i, j;
    bitsliced_gf16_t accumulator;
    bitsliced_addition(&accumulator, &accumulator, &accumulator); //set acc to 0
    for (i = 0; i < 64; i++) {
        for (j = i; j < 64; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            bitsliced_addition(&accumulator, &accumulator, &f.coefficients[position]);
        }
    }
    evaluate_quadratic_polynomials_at_x0_x31(&eval_x0_x31, &f, &x0_x63);
    evaluate_quadratic_polynomials_at_x0_x63(&evaluation, &f, &x0_x63, &eval_x0_x31);
    bitsliced_addition(&evaluation, &evaluation, &accumulator);
    EXPECT_EQ(evaluation.c, 0);
    EXPECT_EQ(evaluation.y, 0);
    EXPECT_EQ(evaluation.x, 0);
    EXPECT_EQ(evaluation.y_x, 0);
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
        prng_gen(&prng, (uint8_t * ) & linear_coefficients, sizeof(bitsliced_gf16_t));
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
    printf("solve_32x32_gf16_system: %lucc\n", total / iterations);
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
    bitsliced_gf16_t x0_x31, x0_x31_prime, y0_y63, evaluation_in_x0_x31;
    bitsliced_quadratic_polynomials_t f;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];

    memset(&x0_x31_prime, 0xff, sizeof(bitsliced_gf16_t));
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t * ) & x0_x31, sizeof(bitsliced_gf16_t));
    x0_x31.c &= 0xFFFFFFFFlu;
    x0_x31.y &= 0xFFFFFFFFlu;
    x0_x31.x &= 0xFFFFFFFFlu;
    x0_x31.y_x &= 0xFFFFFFFFlu;

    int i, j;
    memset(&f, 0x00, sizeof(f));
    for (i = 0; i < 32; i++) {
        for (j = i; j < 64; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            prng_gen(&prng, (uint8_t * ) & f.coefficients[position], sizeof(bitsliced_gf16_t));
        }
    }

    unsigned int dummy;
    unsigned long long t1 = __rdtscp(&dummy);
    for (i = 0; i < 10; i++) {
        find_preimage_of_x0_x31_by_32_polynomials_of_first_layer(&y0_y63, &evaluation_in_x0_x31, &f, &x0_x31, &prng);
    }
    unsigned long long t2 = __rdtscp(&dummy);
    std::cout << "find_preimage_of_x0_x31_by_32_polynomials_of_first_layer: " << (t2 - t1) / 10 << std::endl;
    evaluate_quadratic_polynomials_at_x0_x63(&x0_x31_prime, &f, &y0_y63, &evaluation_in_x0_x31);
    bitsliced_addition(&x0_x31_prime, &x0_x31_prime, &x0_x31);
    EXPECT_EQ(gf16_is_zero(x0_x31_prime, 0), 1);
}

TEST(sign_tests, find_preimage_of_x64_x96_by_32_polynomials_of_second_layer) {
    //find the solution to (f_0(y0, y95), ..., f63(y0, ..., y95)) = (x0, x,63)
    bitsliced_gf16_t x0_x63, x0_x63_prime, y0_y63, y64_y95, evaluation_in_x0_x31;
    bitsliced_quadratic_polynomials_t f;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];

    memset(&x0_x63_prime, 0xff, sizeof(bitsliced_gf16_t));
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t * ) & x0_x63, sizeof(bitsliced_gf16_t));

    int i, j;
    memset(&f, 0x00, sizeof(f));
    //first layer
    for (i = 0; i < 32; i++) {
        for (j = i; j < 64; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            prng_gen(&prng, (uint8_t * ) & f.coefficients[position], sizeof(bitsliced_gf16_t));
        }
        for (j = 64; j < 96; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            prng_gen(&prng, (uint8_t * ) & f.coefficients[position], sizeof(bitsliced_gf16_t));
            shift_left_gf16(&f.coefficients[position], &f.coefficients[position], 32);
        }
    }
    //second layer
    for (i = 32; i < 64; i++) {
        for (j = i; j < 96; j++) {
            int position = i * N - ((i + 1) * i / 2) + j;
            prng_gen(&prng, (uint8_t * ) & f.coefficients[position], sizeof(bitsliced_gf16_t));
            shift_left_gf16(&f.coefficients[position], &f.coefficients[position], 32);
        }
    }

    unsigned int dummy;
    int return_value = 0;
    unsigned long long total = 0;
    for (i = 0; i < 20; i++) {
        unsigned long long t1 = __rdtscp(&dummy);
        find_preimage_of_x0_x31_by_32_polynomials_of_first_layer(&y0_y63, &evaluation_in_x0_x31, &f, &x0_x63, &prng);
        return_value = find_preimage_of_x64_x96_by_32_polynomials_of_second_layer(&y64_y95, &f, &y0_y63, &x0_x63,
                                                                                  &evaluation_in_x0_x31);
        unsigned long long t2 = __rdtscp(&dummy);
        total += t2 - t1;
    }
    printf("find_preimage_of_x64_x96_by_32_polynomials_of_second_layer: %.2fcc\n", total / 20.0);
    evaluate_quadratic_polynomials_of_second_layer_at_x0_x95(&x0_x63_prime, &f, &y0_y63, &y64_y95);
    bitsliced_addition(&x0_x63_prime, &x0_x63_prime, &x0_x63);
    EXPECT_EQ(return_value, 1);
    EXPECT_EQ(gf16_is_zero(x0_x63_prime, 0), 1);
}


TEST(sign_tests, rainbow_sign) {

    private_key_t private_key;

    bitsliced_gf16_t signature[2];
    memset(signature, 0x00, sizeof(signature));

    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    uint8_t message[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    memset(message, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);

    ASSERT_EQ(generate_private_key(&private_key, &prng), SUCCESS);
    EXPECT_EQ(rainbow_sign(signature, &private_key, message, &prng, SECRET_KEY_SEED_BYTE_LENGTH), SUCCESS);
    auto total = time_sign();
    printf("rainbow_sign: %.2fcc\n", total / 20.0);
}

