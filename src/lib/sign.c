//
// Created by Florian Caullery on 8/2/20.
//

#include <memory.h>
#include "sign.h"
#include "keygen.h"
#include "hash_len_config.h"
#include "utils_hash.h"
#include "error_codes.h"

static void print_bitsliced(const bitsliced_gf16_t in, unsigned int bit_position) {
    if ((in[0] >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in[1] >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in[2] >> bit_position) & 0x01u) {
        printf("1");
    } else {
        printf("0");
    }
    if ((in[3] >> bit_position) & 0x01u) {
        printf("1 ");
    } else {
        printf("0 ");
    }
}


void print_32x32_gf16_system(bitsliced_gf16_t m[32]) {
    int i, j;
    for (i = 0; i < 32; i++) {
        printf("l%02d: ", i);
        for (j = 0; j <= 32; j++) {
            print_bitsliced(m[i], j);
        }
        putchar(10);
    }
    putchar(10);
}

#define COND_SWAP(a, b, cond) (((a) ^= ((b) & cond)), ((b) ^= ((a) & cond)), ((a) ^= ((b) & cond)))

int solve_32x32_gf16_system(bitsliced_gf16_t solution, bitsliced_gf16_t equations_coefficients[32]) {

    uint64_t cond = 0, cond1 = 0;
    bitsliced_gf16_t tmp, tmp1, tmp2;
    for (int i = 0; i < 32; i++)
    {
        //constant time line swapping if needed
        cond = 0llu - gf16_is_zero(equations_coefficients[i], i);
        for (int j = i + 1; j < 32; j++)
        {
            cond1 = cond & (~(gf16_is_zero(equations_coefficients[j], i)));
            COND_SWAP(equations_coefficients[i][0], equations_coefficients[j][0], cond1);
            COND_SWAP(equations_coefficients[i][1], equations_coefficients[j][1], cond1);
            COND_SWAP(equations_coefficients[i][2], equations_coefficients[j][2], cond1);
            COND_SWAP(equations_coefficients[i][3], equations_coefficients[j][3], cond1);
            cond ^= cond1;
        }
        if(cond){
            //Did not find a suitable line to swap => all the column is 0
            return SYSTEM_NOT_SOLVABLE;
        }
        expand_bitsliced(tmp, equations_coefficients[i], i);
        bitsliced_inversion(tmp1, tmp);
        bitsliced_multiplication(tmp, tmp1, equations_coefficients[i]); //L_i <- (1/a_i) * L_i
        copy_gf16(equations_coefficients[i], tmp); //L_i <- (1/a_i) * L_i
        for (size_t j = i + 1; j < 32; j++)
        {
            expand_bitsliced(tmp, equations_coefficients[j], i);
            bitsliced_muladd(equations_coefficients[j], tmp, equations_coefficients[i]);
        }
    }

    memset(solution, 0, sizeof(bitsliced_gf16_t));
    for (size_t i = 0; i < 32; i++)
    {
        solution[0] ^= ((equations_coefficients[i][0] >> 32) & 0x01llu) << i;
        solution[1] ^= ((equations_coefficients[i][1] >> 32) & 0x01llu) << i;
        solution[2] ^= ((equations_coefficients[i][2] >> 32) & 0x01llu) << i;
        solution[3] ^= ((equations_coefficients[i][3] >> 32) & 0x01llu) << i;
    }
    for (int i = 31; i >= 0; i--)
    {   
        memset(tmp, 0, sizeof(tmp));
        for (int j = 0; j < i; j++)
        {
            tmp[0] ^= ((equations_coefficients[j][0] >> i) & 0x01llu) << j;
            tmp[1] ^= ((equations_coefficients[j][1] >> i) & 0x01llu) << j;
            tmp[2] ^= ((equations_coefficients[j][2] >> i) & 0x01llu) << j;
            tmp[3] ^= ((equations_coefficients[j][3] >> i) & 0x01llu) << j;
        }
        expand_bitsliced(tmp1, solution, i);
        bitsliced_muladd(solution, tmp1, tmp);
    }
    return SUCCESS;
}


static inline uint64_t bit_parity(uint64_t v){
    v ^= v >> 32;
    v ^= v >> 16;
    v ^= v >> 8;
    v ^= v >> 4;
    v &= 0xf;
    return (0x6996 >> v) & 1;
}

static void prepare_first_layer_polynomial(bitsliced_gf16_t result, first_layer_polynomial_t pol, 
                                            bitsliced_gf16_t x_0__x_v1, bitsliced_gf16_t image, int position){
    
    bitsliced_gf16_t accumulator, acc_tmp, tmp;

    bitsliced_multiplication(acc_tmp, pol[0][0], x_0__x_v1);
    expand_bitsliced(tmp, x_0__x_v1, 0);

    bitsliced_multiplication(accumulator, acc_tmp, tmp);
    bitsliced_multiplication(result, pol[0][1], tmp);

    for (size_t i = 1; i < V1; i++)
    {
        bitsliced_multiplication(acc_tmp, pol[i][0], x_0__x_v1);
        expand_bitsliced(tmp, x_0__x_v1, i);

        bitsliced_muladd(accumulator, acc_tmp, tmp);
        bitsliced_muladd(result, pol[i][1], tmp);
    }

    result[0] ^= (bit_parity(accumulator[0]) << 32) ^ (((image[0] >> position) & 0x01) << 32);
    result[1] ^= (bit_parity(accumulator[1]) << 32) ^ (((image[1] >> position) & 0x01) << 32);
    result[2] ^= (bit_parity(accumulator[2]) << 32) ^ (((image[2] >> position) & 0x01) << 32);
    result[3] ^= (bit_parity(accumulator[3]) << 32) ^ (((image[3] >> position) & 0x01) << 32);
    
}

static void prepare_second_layer_polynomial(bitsliced_gf16_t result, second_layer_polynomial_t pol, 
                                            bitsliced_gf16_t x_0__x_v1, bitsliced_gf16_t x_v1__x_v1_plus_o1,
                                            bitsliced_gf16_t image, int position){

    x_v1__x_v1_plus_o1[0] = (x_v1__x_v1_plus_o1[0] & 0xFFFFFFFFllu) | 0xFFFFFFFF00000000llu;
    x_v1__x_v1_plus_o1[1] = (x_v1__x_v1_plus_o1[1] & 0xFFFFFFFFllu);
    x_v1__x_v1_plus_o1[2] = (x_v1__x_v1_plus_o1[2] & 0xFFFFFFFFllu);
    x_v1__x_v1_plus_o1[3] = (x_v1__x_v1_plus_o1[3] & 0xFFFFFFFFllu);

    bitsliced_gf16_t accumulator[2], acc_tmp[2], tmp;

    bitsliced_multiplication(acc_tmp[0], pol[0][0], x_0__x_v1);
    bitsliced_multiplication(acc_tmp[1], pol[0][1], x_v1__x_v1_plus_o1);
    expand_bitsliced(tmp, x_0__x_v1, 0);

    bitsliced_multiplication(accumulator[0], acc_tmp[0], tmp);
    bitsliced_multiplication(accumulator[1], acc_tmp[1], tmp);

    for (size_t i = 1; i < V1; i++)
    {
        bitsliced_multiplication(acc_tmp[0], pol[i][0], x_0__x_v1);
        bitsliced_multiplication(acc_tmp[1], pol[i][1], x_v1__x_v1_plus_o1);
        expand_bitsliced(tmp, x_0__x_v1, i);

        bitsliced_muladd(accumulator[0], acc_tmp[0], tmp);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {
        bitsliced_multiplication(acc_tmp[1], pol[i][1], x_v1__x_v1_plus_o1);
        expand_bitsliced(tmp, x_v1__x_v1_plus_o1, i - V1);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }

    shift_left_gf16(tmp, accumulator[1], O1);
    bitsliced_addition(tmp, accumulator[0], tmp); 
    shift_right_gf16(result, accumulator[1], O1);
    result[0] ^= (bit_parity(tmp[0]) << 32) ^ (((image[0] >> position) & 0x01) << 32);
    result[1] ^= (bit_parity(tmp[1]) << 32) ^ (((image[1] >> position) & 0x01) << 32);
    result[2] ^= (bit_parity(tmp[2]) << 32) ^ (((image[2] >> position) & 0x01) << 32);
    result[3] ^= (bit_parity(tmp[3]) << 32) ^ (((image[3] >> position) & 0x01) << 32);
    x_v1__x_v1_plus_o1[0] ^= 0xFFFFFFFF00000000llu;
}

void evaluate_central_map(bitsliced_gf16_t image,  central_map_t *f, bitsliced_gf16_t variables[2]){
    bitsliced_gf16_t tmp;
    memset(image, 0, sizeof(bitsliced_gf16_t));
    for (int i = 0; i < O1; i++)
    {
        polynomial_evaluation_first_layer(tmp, f->first_layer_polynomials[i], variables, i);
        bitsliced_addition(image, tmp, image);
    }
    for (int i = 0; i < O2; i++)
    {
        polynomial_evaluation_second_layer(tmp, f->second_layer_polynomials[i], variables, i + O1);
        bitsliced_addition(image, tmp, image);
    }    
}

int find_preimage_by_central_map(bitsliced_gf16_t preimage[2], central_map_t f, bitsliced_gf16_t image){
    
    uint8_t buf[(V1 + 7) / 2];

    int iterations = 100;
    while(iterations != 0)
    {    
        randombytes(buf, sizeof(buf));
        memcpy(&preimage[0][0], buf, (V1 + 7) / 8);
        preimage[0][0] &= (1llu << V1) - 1llu;
        memcpy(&preimage[0][1], buf + (V1 + 7) / 8, (V1 + 7) / 8);
        preimage[0][1] &= (1llu << V1) - 1llu;
        memcpy(&preimage[0][2], buf + (V1 + 7) / 4, (V1 + 7) / 8);
        preimage[0][2] &= (1llu << V1) - 1llu;
        memcpy(&preimage[0][3], buf + ((V1 + 7) / 8) * 3, (V1 + 7) / 8);
        preimage[0][3] &= (1llu << V1) - 1llu;
        
        bitsliced_gf16_t system[32], tmp;
        for (int i = 0; i < O1; i++)
        {
            prepare_first_layer_polynomial(system[i], f.first_layer_polynomials[i], preimage[0], image, i);
        }
        if(solve_32x32_gf16_system(preimage[1], system) != SUCCESS){
            iterations--;
            continue;
        }
        
        for (int i = 0; i < O2; i++)
        {
            prepare_second_layer_polynomial(system[i], f.second_layer_polynomials[i], preimage[0], preimage[1], image, i + O1);
        }
        if(solve_32x32_gf16_system(tmp, system) != SUCCESS){
            iterations--;
            continue;
        }
        shift_left_gf16(tmp, tmp, O1);
        bitsliced_addition(preimage[1], preimage[1], tmp);
        break;
    }
    if(iterations == 0){
        return SYSTEM_NOT_SOLVABLE;
    }else{
        return SUCCESS;
    }
}

int sign_message(bitsliced_gf16_t signature[2], central_map_t *f, 
                bitsliced_gf16_t t_inverse[O1 + V1], bitsliced_gf16_t s[O1],
                uint8_t *message, uint64_t message_length){
    
    uint8_t digest[_HASH_LEN];
    hash_msg(digest, _HASH_LEN, message, message_length);

    bitsliced_gf16_t image, tmp, tmp1[2];
    memcpy(image, digest, sizeof(bitsliced_gf16_t));

    multiply_y_by_s(tmp, s, image);
    if(find_preimage_by_central_map(tmp1, *f, tmp) == SUCCESS){
        multiply_y_by_t(signature, t_inverse, tmp1);
        return SUCCESS;
    }else{
        return SIGNATURE_FAILURE;
    }

}
