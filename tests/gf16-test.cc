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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}