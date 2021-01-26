//
// Created by Florian Caullery on 8/10/20.
//

#include <string.h>
#include "verify.h"
#include "hash_len_config.h"
#include "utils_hash.h"


void evaluate_quadratic_polynomials(bitsliced_gf16_t *evaluations,
                                    bitsliced_quadratic_polynomials_t *f,
                                    bitsliced_gf16_t *x0_x63, bitsliced_gf16_t *x64_x95) {
    bitsliced_gf16_t x0_x95_expanded[N];
    unsigned int i, j;
    for (i = 0; i < 64; i++) {
        x0_x95_expanded[i].c = -1 * ((x0_x63->c >> i) & 0x01u);
        x0_x95_expanded[i].y = -1 * ((x0_x63->y >> i) & 0x01u);
        x0_x95_expanded[i].x = -1 * ((x0_x63->x >> i) & 0x01u);
        x0_x95_expanded[i].y_x = -1 * ((x0_x63->y_x >> i) & 0x01u);
    }
    for (i = 64; i < N; i++) {
        x0_x95_expanded[i].c = -1 * ((x64_x95->c >> (i - 64)) & 0x01u);
        x0_x95_expanded[i].y = -1 * ((x64_x95->y >> (i - 64)) & 0x01u);
        x0_x95_expanded[i].x = -1 * ((x64_x95->x >> (i - 64)) & 0x01u);
        x0_x95_expanded[i].y_x = -1 * ((x64_x95->y_x >> (i - 64)) & 0x01u);
    }
    memset(evaluations, 0x00, sizeof(bitsliced_gf16_t));
    int offset = 0;
    bitsliced_gf16_t tmp, tmp1, tmp2;
    for (i = 0; i < N; i++) {
        bitsliced_square(&tmp1, &x0_x95_expanded[i]);
        bitsliced_multiplication(&tmp2, &tmp1, &f->coefficients[offset]);
        bitsliced_addition(evaluations, evaluations, &tmp2);
        offset++;
        tmp2.c = 0;
        tmp2.y = 0;
        tmp2.x = 0;
        tmp2.y_x = 0;
        for (j = i + 1; j < N; j++) {
            bitsliced_multiplication(&tmp1, &x0_x95_expanded[j], &f->coefficients[offset]);
            bitsliced_addition(&tmp2, &tmp1, &tmp2);
            offset++;
        }
        bitsliced_multiplication(&tmp1, &tmp2, &x0_x95_expanded[i]);
        bitsliced_addition(evaluations, evaluations, &tmp2);
    }
}

int rainbow_verify(bitsliced_gf16_t *signature, public_key_t *public_key, uint8_t *message, uint64_t message_length) {
    int offset = 0;
    bitsliced_gf16_t x, signature_evaluation;
    memset(&x, 0, sizeof(x));
    uint8_t digest[_HASH_LEN];
    hash_msg(digest, _HASH_LEN, message, message_length);
    memcpy(&x.c, digest, sizeof(uint64_t));
    offset += sizeof(uint64_t);
    memcpy(&x.y, digest + offset, sizeof(uint64_t));
    offset += sizeof(uint64_t);
    memcpy(&x.x, digest + offset, sizeof(uint64_t));
    offset += sizeof(uint64_t);
    memcpy(&x.y_x, digest + offset, sizeof(uint64_t));

    evaluate_quadratic_polynomials(&signature_evaluation, &public_key->polynomials, signature, signature + 1);
    bitsliced_addition(&signature_evaluation, &signature_evaluation, &x);
    uint64_t signature_verification =
            signature_evaluation.c | signature_evaluation.y | signature_evaluation.x | signature_evaluation.y_x;
    return signature_verification == 0;
}

