//
// Created by Florian Caullery on 8/10/20.
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
#include "../src/lib/gf16.h"
#include "../src/lib/verify.h"
}



TEST(verify_tests, evaluate_polynomial) {
    bitsliced_gf16_t x0_x63, x64_x95, evaluation;
    bitsliced_quadratic_polynomials_t f;
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];

    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t *) &x0_x63, sizeof(bitsliced_gf16_t));

    int i;
    for (i = 0; i < N * (N + 1) / 2; i++) {
        prng_gen(&prng, (uint8_t *) &f.coefficients[i], sizeof(bitsliced_gf16_t));
    }

    unsigned int dummy;
    unsigned long long t1 = __rdtscp(&dummy);
    for (i = 0; i < 200; i++) {
        evaluate_quadratic_polynomials(&evaluation, &f, &x0_x63, &x64_x95);
    }
    unsigned long long t2 = __rdtscp(&dummy);
    printf("verify: %.2fcc\n", (float) (t2 - t1) / 200.0);
    EXPECT_EQ(gf16_is_zero(x0_x63, 0), 1);
}