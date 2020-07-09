#include "gf16.h"
#include <time.h>
#include <stdio.h>

void bitsliced_addition(bitsliced_gf16_t a_times_b, bitsliced_gf16_t a, bitsliced_gf16_t b) {
    clock_t time = clock();
    printf("%lu", time);
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