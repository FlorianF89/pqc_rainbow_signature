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


uint64_t gf16_is_zero(bitsliced_gf16_t a, unsigned int position);

void copy_gf16(bitsliced_gf16_t *destination, bitsliced_gf16_t *source);

void shift_left_gf16(bitsliced_gf16_t *destination, bitsliced_gf16_t *source, uint8_t shift);

void shift_right_gf16(bitsliced_gf16_t *destination, bitsliced_gf16_t *source, uint8_t shift);

void move_two_halves_gf16_into_one(bitsliced_gf16_t *destination, bitsliced_gf16_t *up_half,
                                   bitsliced_gf16_t *low_half);

void extract_one_gf16_element_and_place_it_in_given_position(bitsliced_gf16_t *destination,
                                                             uint8_t destination_position, bitsliced_gf16_t *in,
                                                             uint8_t in_position);

int bitsliced_gf16_is_one(bitsliced_gf16_t in);

void bitsliced_addition(bitsliced_gf16_t *a_times_b, bitsliced_gf16_t *a, bitsliced_gf16_t *b);

void bitsliced_multiplication(bitsliced_gf16_t *a_times_b, const bitsliced_gf16_t *a, const bitsliced_gf16_t *b);

void
bitsliced_vectorized_multiplication(bitsliced_gf16_t *a_times_b, const bitsliced_gf16_t *a, const bitsliced_gf16_t *b);

void bitsliced_square(bitsliced_gf16_t *a_square, bitsliced_gf16_t *a);

void bitsliced_inversion(bitsliced_gf16_t *a_inverse, bitsliced_gf16_t *a);

#endif
