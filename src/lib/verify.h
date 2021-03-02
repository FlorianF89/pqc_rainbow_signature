//
// Created by Florian Caullery on 8/10/20.
//

#ifndef PQC_RAINBOW_VERIFY_H
#define PQC_RAINBOW_VERIFY_H


#include "gf16.h"
#include "keygen.h"

int rainbow_verify(bitsliced_gf16_t signature[2], 
                    bitsliced_gf16_t public_key[NUMBER_OF_EQUATIONS][NUMBER_OF_VARIABLES][2], 
                    uint8_t *message, uint64_t message_length);

#endif //PQC_RAINBOW_VERIFY_H
