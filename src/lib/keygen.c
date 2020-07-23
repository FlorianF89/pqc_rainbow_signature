//
// Created by Florian Caullery on 7/20/20.
//

#include <memory.h>
#include "keygen.h"
#include "utils_prng.h"
#include "error_codes.h"

void print_binary_matrix(const uint64_t A[64], int is_transposed) {
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

void transpose64(uint64_t *A) {
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
        for (m = 0; m < 64; m++) {
            printf("%llu", (A[j] >> m) & 1u);
        }
        putchar(10);
        A[j] = ((A[j] >> 1u) & 0x5555555555555555llu) | ((A[j] & 0x5555555555555555llu) << 1u);
// swap consecutive pairs
        A[j] = ((A[j] >> 2u) & 0x3333333333333333llu) | ((A[j] & 0x3333333333333333llu) << 2u);
// swap nibbles ...
        A[j] = ((A[j] >> 4u) & 0x0F0F0F0F0F0F0F0Fllu) | ((A[j] & 0x0F0F0F0F0F0F0F0Fllu) << 4u);
// swap bytes
        A[j] = ((A[j] >> 8u) & 0x00FF00FF00FF00FFllu) | ((A[j] & 0x00FF00FF00FF00FFllu) << 8u);
// swap 2-byte long pairs
        A[j] = ((A[j] >> 16u) & 0x0000FFFF0000FFFFllu) | ((A[j] & 0x0000FFFF0000FFFFllu) << 16u);
        A[j] = (A[j] >> 32u) | (A[j] << 32u);
        for (n = 63; n >= 0; n--) {
            printf("%llu", (A[j] >> n) & 1u);
        }
        putchar(10);
        putchar(10);

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
    print_binary_matrix(A, 0);
    transpose64(A);
    print_binary_matrix(A, 1);
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
