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

bitsliced_gf16_t extract_then_expand(bitsliced_gf16_t a, unsigned int position, unsigned int shift) {
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
        for (unsigned int j = i + 1; j < 32; ++j) {
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

uint32_t parity_of_32_bit_words(uint32_t in) {
    in ^= in >> 16u;
    in ^= in >> 8u;
    in ^= in >> 4u;
    in &= 0xfu;
    return (0x6996u >> in) & 1u;
}

uint64_t parity_of_64_bit_words(uint64_t in) {
    in ^= in >> 32u;
    in ^= in >> 16u;
    in ^= in >> 8u;
    in ^= in >> 4u;
    in &= 0xfu;
    return (0x6996u >> in) & 1u;
}

void bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(bitsliced_gf16_t *out,
                                                                             bitsliced_gf16_t *in,
                                                                             uint32_t position_to_place) {

    out->c = (parity_of_32_bit_words(in->c)) << position_to_place;
    out->y = (parity_of_32_bit_words(in->y)) << position_to_place;
    out->x = (parity_of_32_bit_words(in->x)) << position_to_place;
    out->y_x = (parity_of_32_bit_words(in->y_x)) << position_to_place;

}

void bitsliced_gf16_sum_64_first_elements_and_place_result_in_given_position(bitsliced_gf16_t *out,
                                                                             bitsliced_gf16_t *in,
                                                                             uint32_t position_to_place) {

    out->c = (parity_of_64_bit_words(in->c)) << position_to_place;
    out->y = (parity_of_64_bit_words(in->y)) << position_to_place;
    out->x = (parity_of_64_bit_words(in->x)) << position_to_place;
    out->y_x = (parity_of_64_bit_words(in->y_x)) << position_to_place;

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
            bitsliced_gf16_sum_32_first_elements_and_place_result_in_given_position(&tmp, &tmp, i);
            bitsliced_addition(&a_times_b[j], &a_times_b[j], &tmp);
        }
    }
}

void add_32x32_gf16_matrices(bitsliced_gf16_t a_plus_b[32], bitsliced_gf16_t a[32], bitsliced_gf16_t b[32]) {

    uint32_t i;
    for (i = 0; i < 32; i++) {
        bitsliced_addition(&a_plus_b[i], &a[i], &b[i]);
    }
}

int generate_random_32x32_gf16_matrix(bitsliced_gf16_t a[32], prng_t *prng) {
    int i;
    if (prng == NULL) {
        return PRNG_FAILURE;
    }
    for (i = 0; i < 32; i++) {
        prng_gen(prng, (unsigned char *) &a[i].c, sizeof(uint32_t));
        a[i].c &= 0xFFFFFFFFu;
        prng_gen(prng, (unsigned char *) &a[i].y, sizeof(uint32_t));
        a[i].y &= 0xFFFFFFFFu;
        prng_gen(prng, (unsigned char *) &a[i].x, sizeof(uint32_t));
        a[i].x &= 0xFFFFFFFFu;
        prng_gen(prng, (unsigned char *) &a[i].y_x, sizeof(uint32_t));
        a[i].y_x &= 0xFFFFFFFFu;
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

int generate_random_matrix_s(private_key_t *private_key, prng_t *prng) {

//    memset(private_key->s_prime, 0, sizeof(private_key->s_prime));
    if (generate_random_32x32_gf16_matrix(private_key->s_prime, prng) != SUCCESS) {
        return PRNG_FAILURE;
    }

    memset(private_key->s, 0, sizeof(matrix_s_t));

    int i;
    uint64_t j = 1u;
    for (i = 0; i < O1; i++) {
        private_key->s[i].c = j;
        j <<= 1u;
    }
    for (i = O1; i < O1 + O2; i++) {
        copy_gf16(&private_key->s[i], &private_key->s_prime[i - O1]);
        private_key->s[i].c |= j;
        j <<= 1u;
    }
    return SUCCESS;
}


int generate_random_matrix_t(private_key_t *private_key, prng_t *prng) {

    memset(private_key->t1, 0, sizeof(private_key->t1));
    memset(private_key->t2, 0, sizeof(private_key->t2));
    memset(private_key->t3, 0, sizeof(private_key->t3));
    if (generate_random_32x32_gf16_matrix(private_key->t1, prng) != SUCCESS) {
        return PRNG_FAILURE;
    }
    if (generate_random_32x32_gf16_matrix(private_key->t2, prng) != SUCCESS) {
        return PRNG_FAILURE;
    }
    if (generate_random_32x32_gf16_matrix(private_key->t3, prng) != SUCCESS) {
        return PRNG_FAILURE;
    }

    memset(private_key->t, 0, sizeof(matrix_t_t));
    int i;
    uint64_t j = 1u;
    for (i = 0; i < V1; i++) {
        private_key->t[0][i].c = j;
        j <<= 1u;
    }
    for (i = V1; i < V1 + O1; i++) {
        copy_gf16(&private_key->t[0][i], &private_key->t1[i - V1]);
        private_key->t[0][i].c |= j;
        j <<= 1u;
    }

    j = 1u;
    for (i = V1 + O1; i < O1 + O2 + V1; i++) {
        move_two_halves_gf16_into_one(&private_key->t[0][i], &private_key->t2[i - V1 - O1],
                                      &private_key->t3[i - V1 - O1]);
        private_key->t[1][i].c = j;
        j <<= 1u;
    }

    return SUCCESS;
}

static int generate_random_matrices_f_i_for_i_in_v1_v2(private_key_t *private_key, prng_t *prng) {

    uint64_t i, j;
    int return_value = 0;
    for (i = 0; i < O1; i++) {
        return_value += generate_random_32x32_gf16_upper_triangular_matrix(private_key->f1s[i], prng);
        for (j = 0; j < O1; j++) {
            copy_gf16(&private_key->f[i][0][j], &private_key->f1s[i][j]);
        }
        return_value += generate_random_32x32_gf16_matrix(private_key->f2s[i], prng);
        for (j = 0; j < O2; j++) {
            copy_gf16(&private_key->f[i][0][j + O1], &private_key->f2s[i][j]);
        }
    }
    return return_value;
}

static int generate_random_matrices_f_i_for_i_in_v2_n(private_key_t *private_key, prng_t *prng) {

    uint64_t i, j;
    int return_value = 0;
    for (i = O1; i < O1 + O2; i++) {
        return_value += generate_random_32x32_gf16_upper_triangular_matrix(private_key->f1s[i], prng);
        for (j = 0; j < O1; j++) {
            copy_gf16(&private_key->f[i][0][j], &private_key->f1s[i][j]);
        }
        return_value += generate_random_32x32_gf16_matrix(private_key->f2s[i], prng);
        return_value += generate_random_32x32_gf16_upper_triangular_matrix(private_key->f5s[i - O1], prng);
        for (j = 0; j < O2; j++) {
            move_two_halves_gf16_into_one(&private_key->f[i][0][j + O1], &private_key->f2s[i][j],
                                          &private_key->f5s[i - O1][j]);
        }
        return_value += generate_random_32x32_gf16_matrix(private_key->f3s[i - O1], prng);
        return_value += generate_random_32x32_gf16_matrix(private_key->f6s[i - O1], prng);
        for (j = 0; j < O2; j++) {
            move_two_halves_gf16_into_one(&private_key->f[i][0][j + V1 + O1], &private_key->f3s[i - O1][j],
                                          &private_key->f6s[i - O1][j]);
        }
    }
    return return_value;
}


int generate_random_matrices_f(private_key_t *private_key, prng_t *prng) {
    memset(private_key->f, 0, sizeof(matrix_fi_t) * (O1 + O2));

    int return_value = generate_random_matrices_f_i_for_i_in_v1_v2(private_key, prng);
    return_value += generate_random_matrices_f_i_for_i_in_v2_n(private_key, prng);

    return return_value;
}

int generate_private_key(private_key_t *private_key, prng_t *prng) {

    int return_value = generate_random_matrix_s(private_key, prng);
    return_value += generate_random_matrix_t(private_key, prng);
    return_value += generate_random_matrices_f(private_key, prng);

    if (return_value != SUCCESS) {
        return PRNG_FAILURE;
    }

    return SUCCESS;
}

void derive_public_key_from_private_key(public_key_t *public_key, private_key_t *private_key) {

    bitsliced_gf16_t tmp[32];
    bitsliced_gf16_t tmp2[32];
    bitsliced_gf16_t tmp3[32];
    bitsliced_gf16_t tmp4[32];
    bitsliced_gf16_t tmp_single_element;
    bitsliced_gf16_t q2[32];
    bitsliced_gf16_t q3[32];
    bitsliced_gf16_t q5[32];
    bitsliced_gf16_t q6[32];
    bitsliced_gf16_t q9[32];
    bitsliced_gf16_t t1_transpose[32];
    bitsliced_gf16_t t2_transpose[32];
    bitsliced_gf16_t t3_transpose[32];

    transpose_32x32_gf16_matrix(t1_transpose, private_key->t1);
    transpose_32x32_gf16_matrix(t2_transpose, private_key->t2);
    transpose_32x32_gf16_matrix(t3_transpose, private_key->t3);
    memset(public_key->mp, 0x00, sizeof(public_key->mp));
    memset(public_key->mq, 0x00, sizeof(public_key->mq));

    uint64_t i, j, k;
    uint32_t offset;
    for (i = 0; i < O1; i++) {
        offset = 0;
        //Q2 = (F1 + F1^T).T1 + F2
        //tmp = F1 + F1^T
        transpose_and_add_32x32_gf16_matrices(tmp, private_key->f1s[i], private_key->f1s[i]);
        multiply_32x32_gf16_matrices(q2, tmp, private_key->t1);
        add_32x32_gf16_matrices(q2, q2, private_key->f2s[i]);
        //fills the ith row of MQ with the elements of Q1||Q2
        for (j = 0; j < O1; j++) {
            for (k = j; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i,
                                                                        &private_key->f1s[i][k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
            for (k = 0; k < O2; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q2[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
        //Q3 = (F1 +F1^T)·T2 + F2 ·T3
        //tmp2 = F2.T3
        multiply_32x32_gf16_matrices(q3, tmp, private_key->t2);
        multiply_32x32_gf16_matrices(tmp2, private_key->f2s[i], private_key->t3);
        add_32x32_gf16_matrices(q3, q3, tmp2);
        //fills the ith row of MQ with the elements of Q3
        for (j = 0; j < O1; j++) {
            for (k = 0; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q3[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
        //Q9 = UT(T2^T ·F1 ·T2 +T2^T ·F2 ·T3 )
        multiply_32x32_gf16_matrices(tmp3, private_key->f1s[i], private_key->t2);
        multiply_32x32_gf16_matrices(q9, t2_transpose, tmp3);
        multiply_32x32_gf16_matrices(tmp3, t2_transpose, tmp2);
        add_32x32_gf16_matrices(q9, q9, tmp3);
        gaussian_elimination_for_32x32_gf16_matrix(q9);
        //Q6 = T1^T(F1 +F1^T)·T2 +T1^T ·F2 ·T3 +F2^T ·T2
        add_32x32_gf16_matrices(tmp3, tmp, tmp2);
        multiply_32x32_gf16_matrices(q6, t1_transpose, tmp3);
        transpose_32x32_gf16_matrix(tmp3, private_key->f2s[i]);
        multiply_32x32_gf16_matrices(tmp2, tmp3, private_key->t2);
        add_32x32_gf16_matrices(q6, q6, tmp2);
        //Q5 = UT(T1^T · F1 · T1 + T1^T ·F2)
        multiply_32x32_gf16_matrices(tmp, private_key->f1s[i], private_key->t1);
        add_32x32_gf16_matrices(tmp, tmp, private_key->f2s[i]);
        multiply_32x32_gf16_matrices(q5, t1_transpose, tmp);
        gaussian_elimination_for_32x32_gf16_matrix(q5);

        //fills the ith row of MQ with the elements of Q5||Q6
        for (j = 0; j < O1; j++) {
            for (k = j; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q5[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
            for (k = 0; k < O2; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q6[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
        //fills the ith row of MQ with the elements of Q9
        for (j = 0; j < O1; j++) {
            for (k = j; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q9[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
    }
    for (i = O1; i < O1 + O2; i++) {
        offset = 0;
        //Q2 = (F1 + F1^T).T1 + F2
        //tmp = F1 + F1^T
        transpose_and_add_32x32_gf16_matrices(tmp, private_key->f1s[i], private_key->f1s[i]);
        multiply_32x32_gf16_matrices(q2, tmp, private_key->t1);
        add_32x32_gf16_matrices(q2, q2, private_key->f2s[i]);
        //fills the ith row of MQ with the elements of Q1||Q2
        for (j = 0; j < O1; j++) {
            for (k = j; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i,
                                                                        &private_key->f1s[i][k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
            for (k = 0; k < O2; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q2[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
        //Q3 = (F1 +F1^T)·T2 + F2 ·T3 + F3
        //tmp2 = F2.T3
        multiply_32x32_gf16_matrices(q3, tmp, private_key->t2);
        multiply_32x32_gf16_matrices(tmp2, private_key->f2s[i], private_key->t3);
        add_32x32_gf16_matrices(q3, q3, tmp2);
        add_32x32_gf16_matrices(q3, q3, private_key->f3s[i - O1]);
        //fills the ith row of MQ with the elements of Q3
        for (j = 0; j < O1; j++) {
            for (k = 0; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q3[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
        //Q9 = UT(T2^T ·F1 ·T2 +T2^T ·F2 ·T3 + T3^T . F5 . T3 + T2^T . F3 + T3^T . F6)
        multiply_32x32_gf16_matrices(tmp3, private_key->f1s[i], private_key->t2); //tmp3 = F1 T2
        multiply_32x32_gf16_matrices(q9, t2_transpose, tmp3); //q9 = T2^T . F1. T2
        multiply_32x32_gf16_matrices(tmp3, t2_transpose, tmp2); // tmp3 = T2^T F3 T2
        add_32x32_gf16_matrices(q9, q9, tmp3); //q9 = T2^T . F1. T2 + T2^T F3 T2
        multiply_32x32_gf16_matrices(tmp3, private_key->f5s[i - O1], private_key->t3); // tmp3 = F5 T3
        multiply_32x32_gf16_matrices(tmp4, t3_transpose, tmp3); //tmp4 = T3^t F5 T3
        add_32x32_gf16_matrices(q9, q9, tmp4); // q9 = T2^T ·F1 ·T2 +T2^T ·F2 ·T3 + T3^T . F5 . T3
        multiply_32x32_gf16_matrices(tmp3, t2_transpose, private_key->f3s[i - O1]); //tmp3 = T2^T F3
        add_32x32_gf16_matrices(q9, q9, tmp3); // q9 = T2^T ·F1 ·T2 +T2^T ·F2 ·T3 + T3^T . F5 . T3 + T2^T F3
        multiply_32x32_gf16_matrices(tmp3, t3_transpose, private_key->f6s[i - O1]); //tmp3 = T3^T F6
        add_32x32_gf16_matrices(q9, q9, tmp3); // q9 = T2^T ·F1 ·T2 +T2^T ·F2 ·T3 + T3^T . F5 . T3 + T2^T F3 + T3^T . F6
        gaussian_elimination_for_32x32_gf16_matrix(q9);

        //Q6 = T1^T(F1 +F1^T)·T2 +T1^T ·F2 ·T3 +F2^T ·T2 + T1^T F3 + (F5 + F5^T) . T3 + F6
        add_32x32_gf16_matrices(tmp3, tmp, tmp2);
        multiply_32x32_gf16_matrices(q6, t1_transpose, tmp3);
        transpose_32x32_gf16_matrix(tmp3, private_key->f2s[i]);
        multiply_32x32_gf16_matrices(tmp2, tmp3, private_key->t2);
        add_32x32_gf16_matrices(q6, q6, tmp2);
        multiply_32x32_gf16_matrices(tmp3, t1_transpose, private_key->f3s[i - O1]);
        add_32x32_gf16_matrices(q6, q6, tmp3);
        transpose_and_add_32x32_gf16_matrices(tmp3, private_key->f5s[i - O1], private_key->f5s[i - O1]);
        multiply_32x32_gf16_matrices(tmp4, tmp3, private_key->t3);
        add_32x32_gf16_matrices(q6, q6, tmp4);
        add_32x32_gf16_matrices(q6, q6, private_key->f6s[i - O1]);

        //Q5 = UT(T1^T · F1 · T1 + T1^T ·F2 + F5)
        multiply_32x32_gf16_matrices(tmp, private_key->f1s[i], private_key->t1);
        add_32x32_gf16_matrices(tmp, tmp, private_key->f2s[i]);
        multiply_32x32_gf16_matrices(q5, t1_transpose, tmp);
        add_32x32_gf16_matrices(q5, q5, private_key->f5s[i - O1]);
        gaussian_elimination_for_32x32_gf16_matrix(q5);

        //fills the ith row of MQ with the elements of Q5||Q6
        for (j = 0; j < O1; j++) {
            for (k = j; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q5[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
            for (k = 0; k < O2; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q6[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
        //fills the ith row of MQ with the elements of Q9
        for (j = 0; j < O1; j++) {
            for (k = j; k < O1; k++) {
                extract_one_gf16_element_and_place_it_in_given_position(&tmp_single_element, i, &q9[k], j);
                bitsliced_addition(&public_key->mq[offset], &public_key->mq[offset], &tmp_single_element);
                offset++;
            }
        }
    }
    memcpy(public_key->mp, public_key->mq, sizeof(public_key->mq));
    for (i = 0; i < N * (N + 1) / 64; i += 32) {
        for (j = 0; j < 32; j++) {
            copy_gf16(&tmp[j], &public_key->mq[i + j]);
            shift_left_gf16(&tmp[j], &tmp[j], 32);
        }
        multiply_32x32_gf16_matrices(tmp2, private_key->s_prime, tmp);
        add_32x32_gf16_matrices(&public_key->mp[i], &public_key->mp[i], tmp2);
    }
}

static inline int determine_position_of_x_i_times_x_j_in_mq_polynomial(unsigned int i, unsigned int j) {
    if (i < j) {
        return i * N - ((i + 1) * i / 2) + j;
    } else {
        return j * N - ((j + 1) * j / 2) + i;
    }
}

void replace_variable_by_linear_combination_in_quadratic_polynomial(bitsliced_quadratic_polynomials_t *modified_f,
                                                                    bitsliced_quadratic_polynomials_t *original_f,
                                                                    int variable_index,
                                                                    bitsliced_gf16_t *linear_combination) {

    memcpy(modified_f, original_f, sizeof(bitsliced_quadratic_polynomials_t));
    unsigned int i;
    int linear_combination_length, offset;
    bitsliced_gf16_t linear_coefficients[64];
    if (variable_index < 32) {
        for (i = 0; i < 64; i++) {
            linear_coefficients[i].c = ((linear_combination->c >> i) & 0x01u) * (-1);
            linear_coefficients[i].y = ((linear_combination->y >> i) & 0x01u) * (-1);
            linear_coefficients[i].x = ((linear_combination->x >> i) & 0x01u) * (-1);
            linear_coefficients[i].y_x = ((linear_combination->y_x >> i) & 0x01u) * (-1);
        }
        linear_combination_length = 64;
        offset = N - 64;
    } else if (variable_index < 64) {
        for (i = 0; i < 32; i++) {
            linear_coefficients[i].c = ((linear_combination->c >> i) & 0x01u) * (-1);
            linear_coefficients[i].y = ((linear_combination->y >> i) & 0x01u) * (-1);
            linear_coefficients[i].x = ((linear_combination->x >> i) & 0x01u) * (-1);
            linear_coefficients[i].y_x = ((linear_combination->y_x >> i) & 0x01u) * (-1);
        }
        linear_combination_length = 32;
        offset = N - 32;
    } else {
        return;
    }

    i = 0;
    int left_to_replace = N;
    bitsliced_gf16_t tmp, tmp1;
    unsigned int k;
    while (i < variable_index) {
        for (k = 0; k < linear_combination_length; k++) {
            bitsliced_multiplication(&tmp, &linear_coefficients[k],
                                     &original_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i,
                                                                                                                    variable_index)]);
            bitsliced_addition(
                    &modified_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i,
                                                                                                   k + offset)],
                    &modified_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i,
                                                                                                   k + offset)],
                    &tmp);
        }
        i++;
        left_to_replace--;
    }
    for (k = 0; k < linear_combination_length; k++) {
        bitsliced_square(&tmp, &linear_coefficients[k]);
        bitsliced_multiplication(&tmp1, &tmp,
                                 &original_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i, i)]);
        bitsliced_addition(
                &modified_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(k + offset, k + offset)],
                &modified_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(k + offset, k + offset)],
                &tmp);
    }
    left_to_replace--;
    i++;
    while (left_to_replace != 0) {
        for (k = 0; k < linear_combination_length; k++) {
            bitsliced_multiplication(&tmp, &linear_coefficients[k],
                                     &original_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i,
                                                                                                                    variable_index)]);
            bitsliced_addition(
                    &modified_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i, k + offset)],
                    &modified_f->coefficients[determine_position_of_x_i_times_x_j_in_mq_polynomial(i, k + offset)],
                    &tmp);
        }
        i++;
        left_to_replace--;
    }
}

