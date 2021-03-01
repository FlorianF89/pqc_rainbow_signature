//
// Created by Florian Caullery on 8/10/20.
//

#include <string.h>
#include "verify.h"
#include "hash_len_config.h"
#include "utils_hash.h"

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

