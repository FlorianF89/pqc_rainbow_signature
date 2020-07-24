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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}