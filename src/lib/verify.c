//
// Created by Florian Caullery on 8/10/20.
//

#include <string.h>
#include "verify.h"
#include "hash_len_config.h"
#include "utils_hash.h"
#include "error_codes.h"


int rainbow_verify(bitsliced_gf16_t signature[2], 
                    bitsliced_gf16_t public_key[NUMBER_OF_EQUATIONS][NUMBER_OF_VARIABLES][2], 
                    uint8_t *message, uint64_t message_length){

    uint8_t buf[_HASH_LEN];
    bitsliced_gf16_t image, verif;
    
    hash_msg(buf, _HASH_LEN, message, message_length);
    memcpy(image, buf, sizeof(bitsliced_gf16_t));
    int is_valid = 1;
    for (int i = 0; i < NUMBER_OF_EQUATIONS; i++)
    {
        polynomial_evaluation_full(verif, public_key[i], signature, i);
        bitsliced_addition(verif, verif, image);
        is_valid &= gf16_is_zero(verif, i);
    }
    
    return is_valid == 1 ? SUCCESS : SIGNATURE_FAILURE;
}

