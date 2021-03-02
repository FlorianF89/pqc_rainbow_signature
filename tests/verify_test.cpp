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
#include "../src/lib/sign.h"
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
        EXPECT_NO_THROW(evaluate_quadratic_polynomials(&evaluation, &f, &x0_x63, &x64_x95));
    }
    unsigned long long t2 = __rdtscp(&dummy);
    printf("evaluate_quadratic_polynomials: %.2fcc\n", (float) (t2 - t1) / 200.0);
}


TEST(verify_tests, rainbow_verify) {

    private_key_t private_key;
    public_key_t public_key;

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
    EXPECT_NO_THROW(derive_public_key_from_private_key(&public_key, &private_key));

    int i;

    unsigned int dummy;
    int signature_verification;
    unsigned long long total = 0, t1, t2;
    for (i = 0; i < 20; i++) {
        t1 = __rdtscp(&dummy);
        signature_verification = rainbow_verify(signature, &public_key, message, SECRET_KEY_SEED_BYTE_LENGTH);
        t2 = __rdtscp(&dummy);
        total += t2 - t1;
    }
    EXPECT_EQ(signature_verification, 1);
    printf("rainbow_verify: %.2fcc\n", total / 20.0);
}