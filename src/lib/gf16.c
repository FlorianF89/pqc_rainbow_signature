#include "gf16.h"
#include <time.h>
#include <stdio.h>
#include <immintrin.h>

uint64_t bitsliced_gf16_is_one(bitsliced_gf16_t in, int position) {

    uint64_t ret = ((in[0] >> position) & 0x01);
    uint64_t ret2 = ((in[1] >> position) & 0x01) | ((in[2] >> position) & 0x01) | ((in[3] >> position) & 0x01);
    return (uint64_t)((ret == 1) & (ret2 == 0));
}


uint64_t gf16_is_zero(bitsliced_gf16_t a, unsigned int position) {

    uint64_t ret = ((a[0] >> position) & 0x01u) | ((a[1] >> position) & 0x01u) | ((a[2] >> position) & 0x01u) |
                   ((a[3] >> position) & 0x01u);
    return ret == 0;
}

void copy_gf16(bitsliced_gf16_t destination, bitsliced_gf16_t source) {
    destination[0] = source[0];
    destination[1] = source[1];
    destination[2] = source[2];
    destination[3] = source[3];
}

void move_two_halves_gf16_into_one(bitsliced_gf16_t destination, bitsliced_gf16_t up_half,
                                   bitsliced_gf16_t low_half) {
    destination[0] = up_half[0] | (low_half[0] << 32u);
    destination[1] = up_half[1] | (low_half[1] << 32u);
    destination[2] = up_half[2] | (low_half[2] << 32u);
    destination[3] = up_half[3] | (low_half[3] << 32u);
}


void bitsliced_addition(bitsliced_gf16_t a_times_b, bitsliced_gf16_t a, bitsliced_gf16_t b) {
    a_times_b[0] = a[0] ^ b[0];
    a_times_b[1] = a[1] ^ b[1];
    a_times_b[2] = a[2] ^ b[2];
    a_times_b[3] = a[3] ^ b[3];
}

void shift_left_gf16(bitsliced_gf16_t destination, bitsliced_gf16_t source, uint8_t shift) {
    destination[0] = source[0] << shift;
    destination[1] = source[1] << shift;
    destination[2] = source[2] << shift;
    destination[3] = source[3] << shift;
}

void shift_right_gf16(bitsliced_gf16_t destination, bitsliced_gf16_t source, uint8_t shift) {
    destination[0] = source[0] >> shift;
    destination[1] = source[1] >> shift;
    destination[2] = source[2] >> shift;
    destination[3] = source[3] >> shift;
}

void extract_one_gf16_element_and_place_it_in_given_position(bitsliced_gf16_t destination,
                                                             uint8_t destination_position, bitsliced_gf16_t in,
                                                             uint8_t in_position) {
    destination[0] = ((in[0] >> in_position) & 0x01u) << destination_position;
    destination[1] = ((in[1] >> in_position) & 0x01u) << destination_position;
    destination[2] = ((in[2] >> in_position) & 0x01u) << destination_position;
    destination[3] = ((in[3] >> in_position) & 0x01u) << destination_position;
}

inline void bitsliced_multiplication(bitsliced_gf16_t a_times_b, const bitsliced_gf16_t a, const bitsliced_gf16_t b) {
    // a* b = (a[0] & b[0]) ^ (a[1] & b[1]) ^ (a[2] & b[3]) ^ (a[3] & (b[2] ^ b[3])) +
    // e_1 ((a[0] & b[1]) ^ (a[1] & (b[0] ^ b[1])) ^(a[2] & (b[2] ^ b[3])) ^ (a[3] & b[2]))
    a_times_b[0] = (a[0] & b[0]) ^ (a[1] & b[1]) ^ (a[2] & b[3]) ^ (a[3] & (b[2] ^ b[3]));
    a_times_b[1] = (a[0] & b[1]) ^ (a[1] & (b[0] ^ b[1])) ^(a[2] & (b[2] ^ b[3])) ^ (a[3] & b[2]);
    a_times_b[2] = (a[0] & b[2]) ^ (a[1] & b[3]) ^ (a[2] & (b[0] ^ b[2])) ^ (a[3] & (b[1] ^ b[3]));
    a_times_b[3] = (a[0] & b[3]) ^ (a[1] & (b[2] ^ b[3])) ^ (a[2] & (b[1] ^ b[3])) ^ (a[3] & (b[0] ^ b[1] ^ b[2] ^ b[3]));
}

inline void bitsliced_muladd(bitsliced_gf16_t a_times_b_plus_c, const bitsliced_gf16_t a, const bitsliced_gf16_t b) {
    a_times_b_plus_c[0] ^= (a[0] & b[0]) ^ (a[1] & b[1]) ^ (a[2] & b[3]) ^ (a[3] & (b[2] ^ b[3]));
    a_times_b_plus_c[1] ^= (a[0] & b[1]) ^ (a[1] & (b[0] ^ b[1])) ^(a[2] & (b[2] ^ b[3])) ^ (a[3] & b[2]);
    a_times_b_plus_c[2] ^= (a[0] & b[2]) ^ (a[1] & b[3]) ^ (a[2] & (b[0] ^ b[2])) ^ (a[3] & (b[1] ^ b[3]));
    a_times_b_plus_c[3] ^= (a[0] & b[3]) ^ (a[1] & (b[2] ^ b[3])) ^ (a[2] & (b[1] ^ b[3])) ^ (a[3] & (b[0] ^ b[1] ^ b[2] ^ b[3]));
}



static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}


void bitsliced_square(bitsliced_gf16_t a_square, bitsliced_gf16_t a) {

    a_square[0] = a[0] ^ a[3] ^ a[1];
    a_square[1] = a[1] ^ a[2];
    a_square[2] = a[2] ^ a[3];
    a_square[3] = a[3];
}

void bitsliced_inversion(bitsliced_gf16_t a_inverse, bitsliced_gf16_t a) {

    a_inverse[0] = a[0] ^ a[1] ^ a[2];
    a_inverse[1] = a[1] ^ a[2] ^ a[3];
    a_inverse[2] = a[2];
    a_inverse[3] = a[2] ^ a[3];

    a_inverse[0] ^= a[1] & a[3];
    a_inverse[1] ^= a[0] & a[2];

    uint64_t tmp1 = a[0] & a[3];
    a_inverse[0] ^= tmp1;
    a_inverse[2] ^= tmp1;
    a_inverse[3] ^= tmp1;

    uint64_t tmp2 = a[1] & a[2];
    a_inverse[0] ^= tmp2;
    a_inverse[2] ^= tmp2;

    uint64_t tmp3 = tmp1 & a[2];
    a_inverse[0] ^= tmp3;
    a_inverse[2] ^= tmp3;

    tmp3 = tmp1 & a[1];
    a_inverse[0] ^= tmp3;
    a_inverse[1] ^= tmp3;

    tmp3 = tmp2 & a[3];
    a_inverse[0] ^= tmp3;
    a_inverse[1] ^= tmp3;
    a_inverse[2] ^= tmp3;
    a_inverse[3] ^= tmp3;

    a_inverse[0] ^= tmp2 & a[0];
}

inline void bitsliced_conditionnal_addition(bitsliced_gf16_t a_plus_b, 
                                    bitsliced_gf16_t a, 
                                    bitsliced_gf16_t b,
                                    uint64_t cond){

    uint64_t mask = (uint64_t) 0 - (cond & 0x01llu);
    a_plus_b[0] = a[0] ^ (b[0] & mask);
    a_plus_b[1] = a[1] ^ (b[1] & mask);
    a_plus_b[2] = a[2] ^ (b[2] & mask);
    a_plus_b[3] = a[3] ^ (b[3] & mask);
}

inline uint64_t expand_bit(uint64_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu) + 1;
    uint64_t ret =  ((in >> c) & 0x1llu);
    return 0 - ret;
}


inline void expand_bitsliced(bitsliced_gf16_t out, bitsliced_gf16_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu);
    out[0] = 0 - ((in[0] >> c) & 0x1llu);
    out[1] = 0 - ((in[1] >> c) & 0x1llu);
    out[2] = 0 - ((in[2] >> c) & 0x1llu);
    out[3] = 0 - ((in[3] >> c) & 0x1llu);
}

inline void expand_and_add_bitsliced(bitsliced_gf16_t out, bitsliced_gf16_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu);
    out[0] ^= 0 - ((in[0] >> c) & 0x1llu);
    out[1] ^= 0 - ((in[1] >> c) & 0x1llu);
    out[2] ^= 0 - ((in[2] >> c) & 0x1llu);
    out[3] ^= 0 - ((in[3] >> c) & 0x1llu);
}
