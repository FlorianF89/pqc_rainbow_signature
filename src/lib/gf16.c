#include "gf16.h"
#include <time.h>
#include <stdio.h>
#include <immintrin.h>

int bitsliced_gf16_is_one(bitsliced_gf16_t in) {

    if ((in.c == 1u) && (in.y_x == 0) && (in.y == 0) && (in.x == 0)) {
        return 1;
    } else {
        return 0;
    }
}


uint64_t gf16_is_zero(bitsliced_gf16_t a, unsigned int position) {

    uint64_t ret = ((a.c >> position) & 0x01u) | ((a.y >> position) & 0x01u) | ((a.x >> position) & 0x01u) |
                   ((a.y_x >> position) & 0x01u);
    return ret == 0;
}

void copy_gf16(bitsliced_gf16_t *destination, bitsliced_gf16_t *source) {
    destination->c = source->c;
    destination->y = source->y;
    destination->x = source->x;
    destination->y_x = source->y_x;
}

void move_two_halves_gf16_into_one(bitsliced_gf16_t *destination, bitsliced_gf16_t *up_half,
                                   bitsliced_gf16_t *low_half) {
    destination->c = up_half->c | (low_half->c << 32u);
    destination->y = up_half->y | (low_half->y << 32u);
    destination->x = up_half->x | (low_half->x << 32u);
    destination->y_x = up_half->y_x | (low_half->y_x << 32u);
}


void bitsliced_addition(bitsliced_gf16_t *a_times_b, bitsliced_gf16_t *a, bitsliced_gf16_t *b) {
    a_times_b->c = a->c ^ b->c;
    a_times_b->y = a->y ^ b->y;
    a_times_b->x = a->x ^ b->x;
    a_times_b->y_x = a->y_x ^ b->y_x;
}

void shift_left_gf16(bitsliced_gf16_t *destination, bitsliced_gf16_t *source, uint8_t shift) {
    destination->c = source->c << shift;
    destination->y = source->y << shift;
    destination->x = source->x << shift;
    destination->y_x = source->y_x << shift;
}

void shift_right_gf16(bitsliced_gf16_t *destination, bitsliced_gf16_t *source, uint8_t shift) {
    destination->c = source->c >> shift;
    destination->y = source->y >> shift;
    destination->x = source->x >> shift;
    destination->y_x = source->y_x >> shift;
}

void extract_one_gf16_element_and_place_it_in_given_position(bitsliced_gf16_t *destination,
                                                             uint8_t destination_position, bitsliced_gf16_t *in,
                                                             uint8_t in_position) {
    destination->c = ((in->c >> in_position) & 0x01u) << destination_position;
    destination->y = ((in->y >> in_position) & 0x01u) << destination_position;
    destination->x = ((in->x >> in_position) & 0x01u) << destination_position;
    destination->y_x = ((in->y_x >> in_position) & 0x01u) << destination_position;
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

void bitsliced_muladd(bitsliced_gf16_t *a_times_b_plus_c, const bitsliced_gf16_t *a, const bitsliced_gf16_t *b) {
    a_times_b_plus_c->c ^= (a->c & b->c) ^ (a->y & b->y) ^  (a->y_x & b->x) ^ (b->y_x & (a->y_x  ^ a->y ^a->x));
    a_times_b_plus_c->y ^= (a->c & b->y) ^ (a->y ^ (b->c ^ b->y)) ^ (a->x & (b->x ^ b->y_x)) ^ (a->y_x & b->x);
    a_times_b_plus_c->x ^= (a->c & b->x) ^ (a->x & (b->c & b->x)) ^ (a->y_x & (b->y ^ b->y_x));
    a_times_b_plus_c->y_x ^= (a->c & b->y_x) ^ (a->y & (b->x ^ b->y_x)) ^ (a->x & (b->y ^ b->y_x)) ^ (a->y_x & (b->c ^ b->y ^ b->x ^ b->y_x));
}



static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}



void _gf16v_madd_u32(bitsliced_gf16_t *accu_c, const bitsliced_gf16_t *a, uint8_t gf16_b, unsigned _num_byte){

    bitsliced_gf16_t b;

    b.c = 0 - ((gf16_b)& 0x01);
    b.x = 0 - ((gf16_b >> 1)& 0x01);
    b.y = 0 - ((gf16_b >> 2)& 0x01);
    b.y_x = 0 - ((gf16_b >> 3)& 0x01);

    for (int i = 0; i < (_num_byte + 31) / 32; ++i)
    {
        bitsliced_muladd(&accu_c[i], &a[i], &b);
    }
}


void bitsliced_vectorized_multiplication(bitsliced_gf16_t *a_times_b, const bitsliced_gf16_t *a,
                                         const bitsliced_gf16_t *b) {

    __m256i v_a_times_b_c, v_a_times_b_y, v_a_times_b_x, v_a_times_b_y_x;
    __m256i v_a_c = {a[0].c, a[1].c, a[2].c, a[3].c};
    __m256i v_a_y = {a[0].y, a[1].y, a[2].y, a[3].y};
    __m256i v_a_x = {a[0].x, a[1].x, a[2].x, a[3].x};
    __m256i v_a_y_x = {a[0].y_x, a[1].y_x, a[2].y_x, a[3].y_x};
    __m256i v_b_c = {b[0].c, b[1].c, b[2].c, b[3].c};
    __m256i v_b_y = {b[0].y, b[1].y, b[2].y, b[3].y};
    __m256i v_b_x = {b[0].x, b[1].x, b[2].x, b[3].x};
    __m256i v_b_y_x = {b[0].y_x, b[1].y_x, b[2].y_x, b[3].y_x};

//
//    v_a_times_b_c = _mm256_and_si256(v_a_c, v_b_c);
//    v_a_times_b_y = _mm256_xor_si256(_mm256_and_si256(v_a_c, v_b_y), _mm256_and_si256(v_a_y, v_b_c));
//    v_a_times_b_x = _mm256_xor_si256(_mm256_and_si256(v_a_c, v_b_x), _mm256_and_si256(v_a_x, v_b_c));
//    v_a_times_b_y_x = _mm256_xor_si256(
//            _mm256_xor_si256(_mm256_and_si256(v_a_c, v_b_y_x), _mm256_and_si256(v_a_y, v_b_x)),
//            _mm256_xor_si256(_mm256_and_si256(v_a_x, v_b_y), _mm256_and_si256(v_a_y_x, v_b_c)));
//
//    __m256i tmp = _mm256_and_si256(v_a_y, v_b_y);
//    v_a_times_b_c = _mm256_xor_si256(tmp, v_a_times_b_c);
//    v_a_times_b_y = _mm256_xor_si256(tmp, v_a_times_b_y);
//
//    tmp = _mm256_and_si256(v_a_y, v_b_y_x);
//    v_a_times_b_x = _mm256_xor_si256(tmp, v_a_times_b_x);
//    v_a_times_b_y_x = _mm256_xor_si256(tmp, v_a_times_b_y_x);
//
//    tmp = _mm256_and_si256(v_a_x, v_b_x);
//    v_a_times_b_y = _mm256_xor_si256(tmp, v_a_times_b_y);
//    v_a_times_b_x = _mm256_xor_si256(tmp, v_a_times_b_x);
//
//    tmp = _mm256_and_si256(v_a_x, v_b_y_x);
//    v_a_times_b_c = _mm256_xor_si256(tmp, v_a_times_b_c);
//    v_a_times_b_y = _mm256_xor_si256(tmp, v_a_times_b_y);
//    v_a_times_b_y_x = _mm256_xor_si256(tmp, v_a_times_b_y_x);
//
//    tmp = _mm256_and_si256(v_a_y_x, v_b_y);
//    v_a_times_b_x = _mm256_xor_si256(tmp, v_a_times_b_x);
//    v_a_times_b_y_x = _mm256_xor_si256(tmp, v_a_times_b_y_x);
//
//    tmp = _mm256_and_si256(v_a_y_x, v_b_x);
//    v_a_times_b_c = _mm256_xor_si256(tmp, v_a_times_b_c);
//    v_a_times_b_y = _mm256_xor_si256(tmp, v_a_times_b_y);
//    v_a_times_b_y_x = _mm256_xor_si256(tmp, v_a_times_b_y_x);
//
//    tmp = _mm256_and_si256(v_a_y_x, v_b_y_x);
//    v_a_times_b_c = _mm256_xor_si256(tmp, v_a_times_b_c);
//    v_a_times_b_x = _mm256_xor_si256(tmp, v_a_times_b_x);
//    v_a_times_b_y_x = _mm256_xor_si256(tmp, v_a_times_b_y_x);

    a_times_b[0].c = v_a_times_b_c[0];
    a_times_b[1].c = v_a_times_b_c[1];
    a_times_b[2].c = v_a_times_b_c[2];
    a_times_b[3].c = v_a_times_b_c[3];
    a_times_b[0].y = v_a_times_b_y[0];
    a_times_b[1].y = v_a_times_b_y[1];
    a_times_b[2].y = v_a_times_b_y[2];
    a_times_b[3].y = v_a_times_b_y[3];
    a_times_b[0].x = v_a_times_b_x[0];
    a_times_b[1].x = v_a_times_b_x[1];
    a_times_b[2].x = v_a_times_b_x[2];
    a_times_b[3].x = v_a_times_b_x[3];
    a_times_b[0].y_x = v_a_times_b_y_x[0];
    a_times_b[1].y_x = v_a_times_b_y_x[1];
    a_times_b[2].y_x = v_a_times_b_y_x[2];
    a_times_b[3].y_x = v_a_times_b_y_x[3];
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



int main(){

    bitsliced_gf16_t a, c;
    int byte_count = 8;
    FILE *fp;
    fp = fopen("/dev/urandom", "rb");
    fread(&a.c, 1, byte_count, fp);
    fread(&a.x, 1, byte_count, fp);
    fread(&a.y, 1, byte_count, fp);
    fread(&a.y_x, 1, byte_count, fp);
    fclose(fp);
    uint8_t b[4];

	clock_t start_t = clock();
	unsigned long ctr = 0;
	while(((clock() - start_t) / CLOCKS_PER_SEC) < 1){
		_gf16v_madd_u32(&c, &a, b[0], 4);
		_gf16v_madd_u32(&c, &a, b[1], 4);
		_gf16v_madd_u32(&c, &a, b[2], 4);
		_gf16v_madd_u32(&c, &a, b[3], 4);
		ctr++;
	}
	printf("opcount %ld\n", ctr * 256);
    return 0;
}
