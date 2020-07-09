#ifndef LIB_HELLO_TIME_H_
#define LIB_HELLO_TIME_H_

#include <stdint.h>

//gf(16) = (GF(2)[y] / y^2 + y + 1)[x] / (x^2 + x + 1) so an element is represented as c_00 + c_01 * y + c_10 * x + c_11 * x * y
typedef struct bitsliced_gf16 {
    uint64_t c;
    uint64_t y;
    uint64_t x;
    uint64_t y_x;
} bitsliced_gf16_t;

void bitsliced_addition(bitsliced_gf16_t a_times_b, bitsliced_gf16_t a, bitsliced_gf16_t b);

void bitsliced_multiplication(bitsliced_gf16_t *a_times_b, const bitsliced_gf16_t *a, const bitsliced_gf16_t *b);

void bitsliced_square(bitsliced_gf16_t *a_square, bitsliced_gf16_t *a);

#endif
