//
// Created by Florian Caullery on 7/20/20.
//

#include <memory.h>
#include <time.h>
#include "keygen.h"
#include "utils_prng.h"
#include "error_codes.h"

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

void transpose_64x64_binary_matrix(uint64_t *A) {
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

void transpose_32x32_matrix(uint32_t *A) {
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
    transpose_64x64_binary_matrix(expanded_matrix_s2);
    for (i = 0; i < 32; i++) {
        s1_plus_s2t[i].c = (expanded_matrix_s1[i] & 0xFFFFFFFFu) ^ (expanded_matrix_s2[63 - i] & 0xFFFFFFFFu);
        s1_plus_s2t[i].y = (expanded_matrix_s1[i] >> 32u) ^ (expanded_matrix_s2[31 - i] & 0xFFFFFFFFu);
        s1_plus_s2t[i].x = (expanded_matrix_s1[i + 32] & 0xFFFFFFFFu) ^ (expanded_matrix_s2[63 - i] >> 32u);
        s1_plus_s2t[i].y_x = (expanded_matrix_s1[i + 32] >> 32u) ^ (expanded_matrix_s2[31 - i] >> 32u);
    }
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

    clock_t t = clock();
    gaussian_elimination_for_32x32_gf16_matrix(s);
    clock_t work_time = (clock() - t);
    printf("work_time = %lu\n", work_time);
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
    memset(f, 0, sizeof(matrix_fi_t) * (O1 + O2));

    generate_random_matrices_f_i_for_i_in_v1_v2(f, prng);
    generate_random_matrices_f_i_for_i_in_v2_n(f + O1, prng);

    return SUCCESS;
}