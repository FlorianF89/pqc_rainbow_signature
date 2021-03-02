//
// Created by Florian Caullery on 8/2/20.
//

#ifndef PQC_RAINBOW_SIGN_H
#define PQC_RAINBOW_SIGN_H

#include "gf16.h"
#include "utils_prng.h"
#include "keygen.h"


void print_32x32_gf16_system(bitsliced_gf16_t m[32]);

int solve_32x32_gf16_system(bitsliced_gf16_t solution, bitsliced_gf16_t equations_coefficients[32]);

void evaluate_central_map(bitsliced_gf16_t image,  central_map_t *f, bitsliced_gf16_t variables[2]);

int find_preimage_by_central_map(bitsliced_gf16_t preimage[2], central_map_t f, bitsliced_gf16_t image);

int sign_message(bitsliced_gf16_t signature[2], central_map_t *f, 
                bitsliced_gf16_t t_inverse[O1 + V1], bitsliced_gf16_t s[O1],
                uint8_t *message, uint64_t message_length);
                
#endif //PQC_RAINBOW_SIGN_H
