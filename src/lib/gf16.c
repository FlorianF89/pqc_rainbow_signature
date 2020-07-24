#include "gf16.h"
#include <time.h>
#include <stdio.h>

int bitsliced_gf16_is_one(bitsliced_gf16_t in) {

    if ((in.c == 1u) && (in.y_x == 0) && (in.y == 0) && (in.x == 0)) {
        return 1;
    } else {
        return 0;
    }
}


uint64_t gf16_is_zero(bitsliced_gf16_t a, unsigned int position) {

    uint64_t ret = (a.c >> position) | (a.y >> position) | (a.x >> position) | (a.y_x >> position);
    return ret & 0x01u;
}


void bitsliced_addition(bitsliced_gf16_t *a_times_b, bitsliced_gf16_t *a, bitsliced_gf16_t *b) {
    a_times_b->c = a->c ^ b->c;
    a_times_b->y = a->y ^ b->y;
    a_times_b->x = a->x ^ b->x;
    a_times_b->y_x = a->y_x ^ b->y_x;
}

void bitsliced_multiplication(bitsliced_gf16_t *a_times_b, const bitsliced_gf16_t *a, const bitsliced_gf16_t *b) {
    a_times_b->c = (a->c & b->c);
    a_times_b->y = (a->c & b->y) ^ (a->y & b->c);
    a_times_b->x = (a->c & b->x) ^ (a->x & b->c);
    a_times_b->y_x = (a->c & b->y_x) ^ (a->y & b->x) ^ (a->x & b->y) ^ (a->y_x & b->c);

    uint64_t tmp = a->y & b->y;
    a_times_b->c ^= tmp;
    a_times_b->y ^= tmp;

    tmp = a->y & b->y_x;
    a_times_b->x ^= tmp;
    a_times_b->y_x ^= tmp;

    tmp = a->x & b->x;
    a_times_b->y ^= tmp;
    a_times_b->x ^= tmp;

    tmp = a->x & b->y_x;
    a_times_b->c ^= tmp;
    a_times_b->y ^= tmp;
    a_times_b->y_x ^= tmp;

    tmp = a->y_x & b->y;
    a_times_b->x ^= tmp;
    a_times_b->y_x ^= tmp;

    tmp = a->y_x & b->x;
    a_times_b->c ^= tmp;
    a_times_b->y ^= tmp;
    a_times_b->y_x ^= tmp;

    tmp = a->y_x & b->y_x;
    a_times_b->c ^= tmp;
    a_times_b->x ^= tmp;
    a_times_b->y_x ^= tmp;
}

void bitsliced_square(bitsliced_gf16_t *a_square, bitsliced_gf16_t *a) {

    a_square->c = a->c ^ a->y_x ^ a->y;
    a_square->y = a->y ^ a->x;
    a_square->x = a->x ^ a->y_x;
    a_square->y_x = a->y_x;

}

void bitsliced_inversion(bitsliced_gf16_t *a_inverse, bitsliced_gf16_t *a) {

    a_inverse->c = a->c ^ a->y ^ a->x;
    a_inverse->y = a->y ^ a->x ^ a->y_x;
    a_inverse->x = a->x;
    a_inverse->y_x = a->x ^ a->y_x;

    a_inverse->c ^= a->y & a->y_x;
    a_inverse->y ^= a->c & a->x;

    uint64_t tmp1 = a->c & a->y_x;
    a_inverse->c ^= tmp1;
    a_inverse->x ^= tmp1;
    a_inverse->y_x ^= tmp1;

    uint64_t tmp2 = a->y & a->x;
    a_inverse->c ^= tmp2;
    a_inverse->x ^= tmp2;

    uint64_t tmp3 = tmp1 & a->x;
    a_inverse->c ^= tmp3;
    a_inverse->x ^= tmp3;

    tmp3 = tmp1 & a->y;
    a_inverse->c ^= tmp3;
    a_inverse->y ^= tmp3;

    tmp3 = tmp2 & a->y_x;
    a_inverse->c ^= tmp3;
    a_inverse->y ^= tmp3;
    a_inverse->x ^= tmp3;
    a_inverse->y_x ^= tmp3;

    a_inverse->c ^= tmp2 & a->c;
}
