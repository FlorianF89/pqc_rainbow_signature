//
// Created by Florian Caullery on 7/20/20.
//

#include <memory.h>
#include <time.h>
#include "keygen.h"
#include "utils_prng.h"
#include "error_codes.h"


static void print_bitsliced(const bitsliced_gf16_t in, unsigned int bit_position) {
    if ((in.c >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in.y >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in.x >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in.y_x >> bit_position) & 0x01u) {
        printf("1 ");
    } else {
        printf("0 ");
    }
}

void print_32x32_gf16_matrix(bitsliced_gf16_t m[32]) {
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            print_bitsliced(m[j], i);
        }
        putchar(10);
    }
    putchar(10);
}

void print_64x64_binary_matrix(const uint64_t *A, int is_transposed) {
    uint64_t i, j;
    putchar(10);
    if (is_transposed == 0) {
        for (j = 0; j < 64; j++) {
            for (i = 0; i < 64; i++) {
                printf("%llu ", (A[i] >> j) & 1llu);
            }
            putchar(10);
        }
    } else {
        for (j = 0; j < 64; j++) {
            for (i = 0; i < 64; i++) {
                printf("%llu ", (A[63 - i] >> j) & 1llu);
            }
            putchar(10);
        }
    }
}

void transpose_reverse_ordered_64x64_binary_matrix(uint64_t *A) {
    uint64_t j, k;
    uint64_t m, t;
    int n;

    for (j = 32, m = 0x00000000FFFFFFFF; j; j >>= 1, m ^= m << j) {
        for (k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            t = (A[k] ^ (A[k | j] >> j)) & m;
            A[k] ^= t;
            A[k | j] ^= (t << j);
        }
    }
    for (j = 0; j < 64; j++) {
        // swap odd and even bits
        A[j] = ((A[j] >> 1u) & 0x5555555555555555llu) | ((A[j] & 0x5555555555555555llu) << 1u);
        A[j] = ((A[j] >> 2u) & 0x3333333333333333llu) | ((A[j] & 0x3333333333333333llu) << 2u);
        A[j] = ((A[j] >> 4u) & 0x0F0F0F0F0F0F0F0Fllu) | ((A[j] & 0x0F0F0F0F0F0F0F0Fllu) << 4u);
        A[j] = ((A[j] >> 8u) & 0x00FF00FF00FF00FFllu) | ((A[j] & 0x00FF00FF00FF00FFllu) << 8u);
        A[j] = ((A[j] >> 16u) & 0x0000FFFF0000FFFFllu) | ((A[j] & 0x0000FFFF0000FFFFllu) << 16u);
        A[j] = (A[j] >> 32u) | (A[j] << 32u);

    }
}

void transpose_64x64_binary_matrix(uint64_t *out, const uint64_t *in) {
    unsigned int i, j, s, d;

    uint64_t x, y;
    uint64_t masks[6][2] = {
            {0x5555555555555555, 0xAAAAAAAAAAAAAAAA},
            {0x3333333333333333, 0xCCCCCCCCCCCCCCCC},
            {0x0F0F0F0F0F0F0F0F, 0xF0F0F0F0F0F0F0F0},
            {0x00FF00FF00FF00FF, 0xFF00FF00FF00FF00},
            {0x0000FFFF0000FFFF, 0xFFFF0000FFFF0000},
            {0x00000000FFFFFFFF, 0xFFFFFFFF00000000}
    };

    for (i = 0; i < 64; i++)
        out[i] = in[i];

    for (d = 5; d >= 0; d--) {
        s = 1u << d;

        for (i = 0; i < 64; i += s * 2)
            for (j = i; j < i + s; j++) {
                x = (out[j] & masks[d][0]) | ((out[j + s] & masks[d][0]) << s);
                y = ((out[j] & masks[d][1]) >> s) | (out[j + s] & masks[d][1]);

                out[j + 0] = x;
                out[j + s] = y;
            }
    }
}

void transpose_reversed_order_32x32_binary_matrix(uint32_t *A) {
    int j, k;
    unsigned m, t;

    m = 0x0000FFFF;
    for (j = 16; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 32; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k + j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k + j] = A[k + j] ^ (t << j);
        }
    }

    //swap
    for (j = 0; j < 32; j++) {
        A[j] = ((A[j] >> 1u) & 0x55555555u) | ((A[j] & 0x55555555u) << 1u);
        A[j] = ((A[j] >> 2u) & 0x33333333u) | ((A[j] & 0x33333333u) << 2u);
        A[j] = ((A[j] >> 4u) & 0x0F0F0F0Fu) | ((A[j] & 0x0F0F0F0Fu) << 4u);
        A[j] = ((A[j] >> 8u) & 0x00FF00FFu) | ((A[j] & 0x00FF00FFu) << 8u);
        A[j] = (A[j] >> 16u) | (A[j] << 16u);
    }
}

void transpose_32x32_gf16_matrix(bitsliced_gf16_t out[32], bitsliced_gf16_t in[32]) {
    uint64_t expanded_matrix[64];

    int i;
    for (i = 0; i < 32; i++) {
        expanded_matrix[i] = ((in[i].c) & 0xFFFFFFFFu) | ((in[i].y) << 32u);
        expanded_matrix[i + 32] = ((in[i].x) & 0xFFFFFFFFu) | ((in[i].y_x) << 32u);
    }
    transpose_reverse_ordered_64x64_binary_matrix(expanded_matrix);
    for (i = 0; i < 32; i++) {
        out[i].c = expanded_matrix[63 - i] & 0xFFFFFFFFu;
        out[i].y = expanded_matrix[31 - i] & 0xFFFFFFFFu;
        out[i].x = expanded_matrix[63 - i] >> 32u;
        out[i].y_x = expanded_matrix[31 - i] >> 32u;
    }
}

static bitsliced_gf16_t extract_then_expand(bitsliced_gf16_t a, unsigned int position, unsigned int shift) {
    bitsliced_gf16_t ret;
    ret.c = ((a.c >> position) & 0x01u) * (0xFFFFFFFFFFFFFFFF << shift);
    ret.x = ((a.x >> position) & 0x01u) * (0xFFFFFFFFFFFFFFFF << shift);
    ret.y = ((a.y >> position) & 0x01u) * (0xFFFFFFFFFFFFFFFF << shift);
    ret.y_x = ((a.y_x >> position) & 0x01u) * (0xFFFFFFFFFFFFFFFF << shift);
    return ret;
}

void gaussian_elimination_for_32x32_gf16_matrix(matrix_s_t s) {
    for (unsigned int i = 0; i < 31; i++) {
        bitsliced_gf16_t tmp = extract_then_expand(s[i], i, i);
        bitsliced_gf16_t tmp1;
        bitsliced_inversion(&tmp1, &tmp);
        bitsliced_multiplication(&tmp, &tmp1, &s[i]);
        for (unsigned int j = i + 1u; j < 32; ++j) {
            bitsliced_gf16_t tmp2 = extract_then_expand(s[j], i, i + 1);
            bitsliced_multiplication(&tmp1, &tmp, &tmp2);
            bitsliced_addition(&s[j], &tmp1, &s[j]);
        }
        s[i].c &= (1u << (i + 1u)) - 1u;
        s[i].y &= (1u << (i + 1u)) - 1u;
        s[i].x &= (1u << (i + 1u)) - 1u;
        s[i].y_x &= (1u << (i + 1u)) - 1u;
    }
}

void transpose_and_add_32x32_gf16_matrices(matrix_s_t s1_plus_s2t, matrix_s_t s1, matrix_s_t s2) {
    uint64_t expanded_matrix_s1[64];
    uint64_t expanded_matrix_s2[64];

    int i;
    for (i = 0; i < 32; i++) {
        expanded_matrix_s1[i] = ((s1[i].c) & 0xFFFFFFFFu) | ((s1[i].y) << 32u);
        expanded_matrix_s2[i] = ((s2[i].c) & 0xFFFFFFFFu) | ((s2[i].y) << 32u);
        expanded_matrix_s1[i + 32] = ((s1[i].x) & 0xFFFFFFFFu) | ((s1[i].y_x) << 32u);
        expanded_matrix_s2[i + 32] = ((s2[i].x) & 0xFFFFFFFFu) | ((s2[i].y_x) << 32u);
    }
    transpose_reverse_ordered_64x64_binary_matrix(expanded_matrix_s2);
    for (i = 0; i < 32; i++) {
        s1_plus_s2t[i].c = (expanded_matrix_s1[i] & 0xFFFFFFFFu) ^ (expanded_matrix_s2[63 - i] & 0xFFFFFFFFu);
        s1_plus_s2t[i].y = (expanded_matrix_s1[i] >> 32u) ^ (expanded_matrix_s2[31 - i] & 0xFFFFFFFFu);
        s1_plus_s2t[i].x = (expanded_matrix_s1[i + 32] & 0xFFFFFFFFu) ^ (expanded_matrix_s2[63 - i] >> 32u);
        s1_plus_s2t[i].y_x = (expanded_matrix_s1[i + 32] >> 32u) ^ (expanded_matrix_s2[31 - i] >> 32u);
    }
}

static uint32_t parity_of_32_bit_words(uint32_t in) {
    in ^= in >> 16u;
    in ^= in >> 8u;
    in ^= in >> 4u;
    in &= 0xfu;
    return (0x6996u >> in) & 1u;
}

static void bitsliced_gf_16_sum_32_first_elements_and_place_result_in_given_position(bitsliced_gf16_t *out,
                                                                                     bitsliced_gf16_t *in,
                                                                                     uint32_t position_to_place) {

    out->c = (parity_of_32_bit_words(in->c)) << position_to_place;
    out->y = (parity_of_32_bit_words(in->y)) << position_to_place;
    out->x = (parity_of_32_bit_words(in->x)) << position_to_place;
    out->y_x = (parity_of_32_bit_words(in->y_x)) << position_to_place;

}

void set_32x32_gf16_matrix_to_identity(bitsliced_gf16_t a[32]) {

    memset(a, 0, 32 * sizeof(bitsliced_gf16_t));
    for (unsigned int i = 0; i < 32; i++) {
        a[i].c = 1u << i;
    }
}

void multiply_32x32_gf16_matrices(bitsliced_gf16_t a_times_b[32], bitsliced_gf16_t a[32], bitsliced_gf16_t b[32]) {

    bitsliced_gf16_t a_transposed[32];
    bitsliced_gf16_t tmp;
    transpose_32x32_gf16_matrix(a_transposed, a);
    memset(a_times_b, 0, 32 * sizeof(bitsliced_gf16_t));
    uint32_t i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            bitsliced_multiplication(&tmp, &a_transposed[i], &b[j]);
            bitsliced_gf_16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, i);
            bitsliced_addition(&a_times_b[j], &a_times_b[j], &tmp);
        }
    }
}

int generate_random_32x32_gf16_matrix(bitsliced_gf16_t a[32], prng_t *prng) {
    int i;
    if (prng == NULL) {
        return PRNG_FAILURE;
    }
    for (i = 0; i < 32; i++) {
        prng_gen(prng, (unsigned char *) &a[i].c, sizeof(uint32_t));
        prng_gen(prng, (unsigned char *) &a[i].y, sizeof(uint32_t));
        prng_gen(prng, (unsigned char *) &a[i].x, sizeof(uint32_t));
        prng_gen(prng, (unsigned char *) &a[i].y_x, sizeof(uint32_t));
    }
    return SUCCESS;
}

int generate_random_32x32_gf16_upper_triangular_matrix(bitsliced_gf16_t a[32], prng_t *prng) {
    unsigned int i;
    if (prng == NULL) {
        return PRNG_FAILURE;
    }
    uint64_t mask = 1u;
    for (i = 0; i < 32; i++) {
        prng_gen(prng, (unsigned char *) &a[i].c, (i + 7) / 8);
        prng_gen(prng, (unsigned char *) &a[i].y, (i + 7) / 8);
        prng_gen(prng, (unsigned char *) &a[i].x, (i + 7) / 8);
        prng_gen(prng, (unsigned char *) &a[i].y_x, (i + 7) / 8);
        a[i].c &= mask;
        a[i].y &= mask;
        a[i].x &= mask;
        a[i].y_x &= mask;
        mask <<= 1u;
        mask |= 1u;
    }
    return SUCCESS;
}

int generate_random_matrix_s(matrix_s_t s, prng_t *prng) {
    if (prng == NULL) {
        return PRNG_FAILURE;
    }
    memset(s, 0, sizeof(matrix_s_t));
    int i;
    uint64_t j = 1u;
    for (i = 0; i < O1 + O2; i++) {
        s[i].c = j;
        j <<= 1u;
    }

    for (i = 0; i < O2; i++) {
        prng_gen(prng, (unsigned char *) &s[O1 + i].c, O1 / 8);
        prng_gen(prng, (unsigned char *) &s[O1 + i].y, O1 / 8);
        prng_gen(prng, (unsigned char *) &s[O1 + i].x, O1 / 8);
        prng_gen(prng, (unsigned char *) &s[O1 + i].y_x, O1 / 8);
    }

    gaussian_elimination_for_32x32_gf16_matrix(s);
    return SUCCESS;
}


int generate_random_matrix_t(matrix_t_t t, prng_t *prng) {

    if (prng == NULL) {
        return PRNG_FAILURE;
    }
    memset(t, 0, sizeof(matrix_t_t));
    int i;
    uint64_t j = 1u;
    for (i = 0; i < V1 + O1; i++) {
        t[0][i].c = j;
        j <<= 1u;
    }
    j = 1u;
    for (i = V1 + O1; i < O1 + O2 + V1; i++) {
        t[1][i].c = j;
        j <<= 1u;
    }

    for (i = 0; i < O1; i++) {
        prng_gen(prng, (unsigned char *) &t[0][V1 + i].c, O1 / 8);
        prng_gen(prng, (unsigned char *) &t[0][V1 + i].y, O1 / 8);
        prng_gen(prng, (unsigned char *) &t[0][V1 + i].x, O1 / 8);
        prng_gen(prng, (unsigned char *) &t[0][V1 + i].y_x, O1 / 8);
    }

    for (i = 0; i < O2; i++) {
        prng_gen(prng, (unsigned char *) &t[0][V1 + O1 + i].c, (V1 + O1) / 8);
        prng_gen(prng, (unsigned char *) &t[0][V1 + O1 + i].y, (V1 + O1) / 8);
        prng_gen(prng, (unsigned char *) &t[0][V1 + O1 + i].x, (V1 + O1) / 8);
        prng_gen(prng, (unsigned char *) &t[0][V1 + O1 + i].y_x, (V1 + O1) / 8);
    }
    return SUCCESS;
}

static int generate_random_matrices_f_i_for_i_in_v1_v2(matrix_fi_t f[O1], prng_t *prng) {

    uint64_t i, j;
    for (i = 0; i < O1; i++) {
        for (j = 0; j < O1 + O2; j++) {
            prng_gen(prng, (unsigned char *) &f[i][0][j].c, O1 / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].y, O1 / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].x, O1 / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].y_x, O1 / 8);
        }
        uint64_t mask = 0x01u;
        for (j = 0; j < O1; j++) {
            f[i][0][j].c &= mask;
            f[i][0][j].y &= mask;
            f[i][0][j].x &= mask;
            f[i][0][j].y_x &= mask;
            mask = (mask << 1u) | 1u;
        }
    }
    uint64_t A[64];
    for (i = 0; i < 64; i++) {
        A[i] = f[0][0][i].c;
    }
    return SUCCESS;
}

static int generate_random_matrices_f_i_for_i_in_v2_n(matrix_fi_t f[O2], prng_t *prng) {

    uint64_t i, j;
    for (i = 0; i < O2; i++) {
        for (j = 0; j < O1 + O2; j++) {
            prng_gen(prng, (unsigned char *) &f[i][0][j].c, O1 / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].y, O1 / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].x, O1 / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].y_x, O1 / 8);
        }
        uint64_t mask = 0x01lu;
        for (j = 0; j < O1 + O2; j++) {
            f[i][0][j].c &= mask;
            f[i][0][j].y &= mask;
            f[i][0][j].x &= mask;
            f[i][0][j].y_x &= mask;
            mask = (mask << 1u) | 1u;
        }
        for (j = O1 + O2; j < O1 + O2 + V1; j++) {
            prng_gen(prng, (unsigned char *) &f[i][0][j].c, (O1 + O2) / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].y, (O1 + O2) / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].x, (O1 + O2) / 8);
            prng_gen(prng, (unsigned char *) &f[i][0][j].y_x, (O1 + O2) / 8);
        }
    }
    return SUCCESS;
}


int generate_random_matrices_f(matrix_fi_t f[O1 + O2], prng_t *prng) {
    if (prng == NULL) {
        return PRNG_FAILURE;
    }
    memset(f, 0, sizeof(matrix_fi_t) * (O1 + O2));

    generate_random_matrices_f_i_for_i_in_v1_v2(f, prng);
    generate_random_matrices_f_i_for_i_in_v2_n(f + O1, prng);

    return SUCCESS;
}

int generate_private_key(private_key_t *private_key, prng_t *prng) {
    if (prng == NULL) {
        return PRNG_FAILURE;
    }

    generate_random_matrix_s(private_key->s, prng);
    generate_random_matrix_t(private_key->t, prng);
    generate_random_matrices_f(private_key->f, prng);

    return SUCCESS;
}
