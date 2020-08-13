#include <x86intrin.h>
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

TEST(gf16_mul, gf16_multiplication_against_test_vectors) {
    bitsliced_gf16_t i, j, i_times_j, result;

    i.c = 0xFFFF0000FFFF0000u;
    i.y = 0xFFFFFFFF00000000u;
    i.x = 0;
    i.y_x = 0;
    j.c = 0xAAAAAAAAAAAAAAAAu;
    j.y = 0xCCCCCCCCCCCCCCCCu;
    j.x = 0xf0f0f0f0f0f0f0f0u;
    j.y_x = 0xff00ff00ff00ff00u;

    i_times_j.c = 0x6666CCCCAAAA0000u;
    i_times_j.y = 0xAAAA6666CCCC0000u;
    i_times_j.x = 0xff0ff00f0f00000u;
    i_times_j.y_x = 0xf0f00ff0ff000000u;

    bitsliced_multiplication(&result, &i, &j);

    EXPECT_EQ (result.c, i_times_j.c);
    EXPECT_EQ (result.y, i_times_j.y);
    EXPECT_EQ (result.x, i_times_j.x);
    EXPECT_EQ (result.y_x, i_times_j.y_x);
}

TEST(gf16_square, gf16_square_against_test_vectors) {
    bitsliced_gf16_t i, i_square, result;

    i.c = 0xAAAA;
    i.y = 0xCCCC;
    i.x = 0xF0F0;
    i.y_x = 0xFF00;

    i_square.c = 0x9966;
    i_square.y = 0x3C3C;
    i_square.x = 0xFF0;
    i_square.y_x = 0xFF00;

    bitsliced_square(&result, &i);

    EXPECT_EQ (result.c, i_square.c);
    EXPECT_EQ (result.y, i_square.y);
    EXPECT_EQ (result.x, i_square.x);
    EXPECT_EQ (result.y_x, i_square.y_x);
}

TEST(gf16_inversion, gf16_inversion_extensive) {
    bitsliced_gf16_t i, i_inverse, ones, result;

    //test all elements of GF16
    i.c = 0xAAAA;
    i.y = 0xCCCC;
    i.x = 0xF0F0;
    i.y_x = 0xFF00;

    ones.c = 0xFFFE;
    ones.y = 0;
    ones.x = 0;
    ones.y_x = 0;

    bitsliced_inversion(&i_inverse, &i);
    bitsliced_multiplication(&result, &i, &i_inverse);

    EXPECT_EQ (result.c, ones.c);
    EXPECT_EQ (result.y, ones.y);
    EXPECT_EQ (result.x, ones.x);
    EXPECT_EQ (result.y_x, ones.y_x);
}

TEST(gf16, compare_vectorized) {
    bitsliced_gf16_t a[4], b[4], a_times_b[4];

    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t *) &a, sizeof(a));
    prng_gen(&prng, (uint8_t *) &b, sizeof(b));
    unsigned int dummy, i;
    unsigned long long total = 0, t1, t2;
    for (i = 0; i < 1000000; i++) {
        t1 = __rdtscp(&dummy);
        bitsliced_multiplication(&a_times_b[0], &a[0], &b[0]);
        bitsliced_multiplication(&a_times_b[1], &a[1], &b[1]);
        bitsliced_multiplication(&a_times_b[2], &a[2], &b[2]);
        bitsliced_multiplication(&a_times_b[3], &a[3], &b[3]);
        t2 = __rdtscp(&dummy);
        total += t2 - t1;
    }
    printf("non vec mul: %.2fcc\n", total / 1000000.0);
    total = 0;
    for (i = 0; i < 1000000; i++) {
        t1 = __rdtscp(&dummy);
        bitsliced_vectorized_multiplication(a_times_b, a, b);
        t2 = __rdtscp(&dummy);
        total += t2 - t1;
    }
    printf("vec mul: %.2fcc\n", total / 1000000.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}