//
// Created by Florian Caullery on 7/20/20.
//

#include <memory.h>
#include <time.h>
#include "keygen.h"
#include "rng.h"
#include "utils_prng.h"
#include "error_codes.h"

#define ELMTS_IN_BITSLICED 64

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


inline uint64_t expand_bit(uint64_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu) + 1;
    uint64_t ret =  ((in >> c) & 0x1llu);
    return 0 - ret;
}


static inline void expand_bitsliced(bitsliced_gf16_t out, bitsliced_gf16_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu);
    out[0] = 0 - ((in[0] >> c) & 0x1llu);
    out[1] = 0 - ((in[1] >> c) & 0x1llu);
    out[2] = 0 - ((in[2] >> c) & 0x1llu);
    out[3] = 0 - ((in[3] >> c) & 0x1llu);
}

static inline void expand_and_add_bitsliced(bitsliced_gf16_t out, bitsliced_gf16_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu);
    out[0] ^= 0 - ((in[0] >> c) & 0x1llu);
    out[1] ^= 0 - ((in[1] >> c) & 0x1llu);
    out[2] ^= 0 - ((in[2] >> c) & 0x1llu);
    out[3] ^= 0 - ((in[3] >> c) & 0x1llu);
}

static inline unsigned int get_square_index_in_polynomial(int variable_index, int number_of_variables){
    return variable_index * number_of_variables - ((variable_index * (variable_index - 1)) / 2);
}


#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
// Assumes i =< j!
static inline unsigned int get_binomial_index_in_polynomial(int i, int j, int number_of_variables){
    if(i > j) SWAP(i, j);
    return i *number_of_variables - ((i * (i - 1)) / 2) + (j - i);
}

void variable_substitution(bitsliced_gf16_t transformed[6], bitsliced_gf16_t original[6], bitsliced_gf16_t variable_list[3]){
    int i, j, k, l;
    bitsliced_gf16_t tmp, tmp2;
    memset(transformed, 0, sizeof(bitsliced_gf16_t) * 6);

    for (i = 0; i < 3; i++)
    {
        bitsliced_square(tmp2, variable_list[i]);
        for(j = i; j < 3; j++){
            expand_bitsliced(tmp, tmp2, j);
            bitsliced_muladd(transformed[get_square_index_in_polynomial(j, 3)], original[get_square_index_in_polynomial(i, 3)], tmp);
        }
        for (j = i + 1; j < 3; j++)
        {
            for(k = i; k < 3; k++){
                expand_bitsliced(tmp2, variable_list[i], k);
                bitsliced_multiplication(tmp, tmp2, original[get_binomial_index_in_polynomial(i, j, 3)]);
                for (l = j; l < 3; l++)
                {
                    expand_bitsliced(tmp2, variable_list[j], l);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(k, l, 3)], tmp, tmp2);
                }
            }
        }
    }
}

#define NUMB_VARS_V1_CHANGE (((((O1 + O2)*(O1 + O2 + 1)) / 2) + ELMTS_IN_BITSLICED - 1) / ELMTS_IN_BITSLICED)

static void inline bitsliced_ror(bitsliced_gf16_t in, unsigned int r){
    uint8_t c = r & 0x3Fu;
    in[0] = (in[0] >> c) | (in[0] << (unsigned) (ELMTS_IN_BITSLICED - c));
    in[1] = (in[1] >> c) | (in[1] << (unsigned) (ELMTS_IN_BITSLICED - c));
    in[2] = (in[2] >> c) | (in[2] << (unsigned) (ELMTS_IN_BITSLICED - c));
    in[3] = (in[3] >> c) | (in[3] << (unsigned) (ELMTS_IN_BITSLICED - c));
}

static void inline bitsliced_rol(bitsliced_gf16_t in, unsigned int r){
    uint8_t c = r & 0x3Fu;
    in[0] = (in[0] << c) | (in[0] >> (unsigned) (ELMTS_IN_BITSLICED - c));
    in[1] = (in[1] << c) | (in[1] >> (unsigned) (ELMTS_IN_BITSLICED - c));
    in[2] = (in[2] << c) | (in[2] >> (unsigned) (ELMTS_IN_BITSLICED - c));
    in[3] = (in[3] << c) | (in[3] >> (unsigned) (ELMTS_IN_BITSLICED - c));
}

static void inline bitsliced_rol_32(bitsliced_gf16_t in, unsigned int r){
    uint8_t c = r & 0x3Fu;

    in[0] = ((in[0] & 0xFFFFFFFF00000000u) << c) | ((in[0] >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c)) & 0xFFFFFFFF00000000u)
            | ((in[0] << c) & 0xFFFFFFFFu) | ((in[0] & 0xFFFFFFFFu) >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c));
    in[1] = ((in[1] & 0xFFFFFFFF00000000u) << c) | ((in[1] >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c)) & 0xFFFFFFFF00000000u)
            | ((in[1] << c) & 0xFFFFFFFFu) | ((in[1] & 0xFFFFFFFFu) >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c));
    in[2] = ((in[2] & 0xFFFFFFFF00000000u) << c) | ((in[2] >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c)) & 0xFFFFFFFF00000000u)
            | ((in[2] << c) & 0xFFFFFFFFu) | ((in[2] & 0xFFFFFFFFu) >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c));
    in[3] = ((in[3] & 0xFFFFFFFF00000000u) << c) | ((in[3] >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c)) & 0xFFFFFFFF00000000u)
            | ((in[3] << c) & 0xFFFFFFFFu) | ((in[3] & 0xFFFFFFFFu) >> (unsigned) ((ELMTS_IN_BITSLICED / 2) - c));
}

//Multiply and put in place for variale change for x_i x_j with i,j in V1
// Writing x_i -> x_i + a_o1 x_o1 + ... + a_{n - 1} x_{n - 1} and
// x_j -> x_j + b_o1 x_o1 + ... + b_{n - 1} x_{n - 1} 
// c will contain in order the coeffcients of [x_o1^2, ... , x_{n - 1}^2, 
// x_o1 x*{o1 + 1}, x_{o1+1} x*{o1+2}, ... , x_{n-1} x*{o1}, x_o1 x*{o1 + 2}, x_{o1+1} x*{o1+3}, ... , x_{n-1} x*{o1 + 1}, ...]
// (rotating fashion)
static void inline gf16_multiply_for_variable_change_v1_v1(bitsliced_gf16_t c[NUMB_VARS_V1_CHANGE], 
                                                            bitsliced_gf16_t a, bitsliced_gf16_t b){
    
    bitsliced_multiplication(c[0], a, b);
    for (unsigned int var_change_idx = 1; var_change_idx < (ELMTS_IN_BITSLICED / 2) ; var_change_idx++)
    {
        bitsliced_ror(b, 1);
        bitsliced_multiplication(c[var_change_idx], a, b);
    }
    //Treats the special case of the rotation by (ELMTS_IN_BITSLICED / 2)
    bitsliced_ror(b, 1);
    bitsliced_multiplication(c[(ELMTS_IN_BITSLICED / 2)], a, b);
    c[(ELMTS_IN_BITSLICED / 2)][0] ^= c[(ELMTS_IN_BITSLICED / 2)][0] >> (ELMTS_IN_BITSLICED / 2);
    c[(ELMTS_IN_BITSLICED / 2)][1] ^= c[(ELMTS_IN_BITSLICED / 2)][1] >> (ELMTS_IN_BITSLICED / 2);
    c[(ELMTS_IN_BITSLICED / 2)][2] ^= c[(ELMTS_IN_BITSLICED / 2)][2] >> (ELMTS_IN_BITSLICED / 2);
    c[(ELMTS_IN_BITSLICED / 2)][3] ^= c[(ELMTS_IN_BITSLICED / 2)][3] >> (ELMTS_IN_BITSLICED / 2);
    
    // In this loop we have to rotate back tmp to align the coefficients with the first loop
    bitsliced_gf16_t tmp_v1_v1;
    for (unsigned int var_change_idx = (ELMTS_IN_BITSLICED / 2) - 1; var_change_idx > 0; var_change_idx--)
    {
        bitsliced_ror(b, 1);
        bitsliced_multiplication(tmp_v1_v1, a, b);
        bitsliced_rol(tmp_v1_v1, var_change_idx);
        bitsliced_addition(c[var_change_idx], c[var_change_idx], tmp_v1_v1);
    }
    //Restores b to its original value
    bitsliced_ror(b, 1);
}

#define NUMB_VARS_O1_CHANGE (((((1 + O2)*(O2))) + ELMTS_IN_BITSLICED - 1) / ELMTS_IN_BITSLICED)

static void inline gf16_multiply_for_variable_change_v1_o1(bitsliced_gf16_t c[NUMB_VARS_V1_CHANGE], 
                                                            bitsliced_gf16_t a, bitsliced_gf16_t b){
    
    bitsliced_gf16_t tmp_v1_o1;
    tmp_v1_o1[0] = (b[0] >> 32u) | (b[0] & 0xFFFFFFFF00000000u);
    tmp_v1_o1[1] = (b[1] >> 32u) | (b[1] & 0xFFFFFFFF00000000u);
    tmp_v1_o1[2] = (b[2] >> 32u) | (b[2] & 0xFFFFFFFF00000000u);
    tmp_v1_o1[3] = (b[3] >> 32u) | (b[3] & 0xFFFFFFFF00000000u);

    bitsliced_multiplication(c[0], a, tmp_v1_o1);
    for (unsigned int var_change_idx = 1; var_change_idx < (ELMTS_IN_BITSLICED / 2) ; var_change_idx++)
    {
        bitsliced_ror(tmp_v1_o1, 1);
        bitsliced_multiplication(c[var_change_idx], a, tmp_v1_o1);
    }
}


static void inline gf16_double_multiply_for_variable_change_o1_o1(bitsliced_gf16_t c[NUMB_VARS_O1_CHANGE], 
                                                            bitsliced_gf16_t a, bitsliced_gf16_t b){
    
    bitsliced_gf16_t tmp_o1_o1;
    tmp_o1_o1[0] = (b[0] >> 32u) | (b[0] & 0xFFFFFFFF00000000u);
    tmp_o1_o1[1] = (b[1] >> 32u) | (b[1] & 0xFFFFFFFF00000000u);
    tmp_o1_o1[2] = (b[2] >> 32u) | (b[2] & 0xFFFFFFFF00000000u);
    tmp_o1_o1[3] = (b[3] >> 32u) | (b[3] & 0xFFFFFFFF00000000u);
    
    unsigned int var_change_idx;
    bitsliced_multiplication(c[0], a, tmp_o1_o1);
    for (var_change_idx = 1; var_change_idx < (ELMTS_IN_BITSLICED / 4) ; var_change_idx++)
    {
        bitsliced_ror(tmp_o1_o1, 1);
        bitsliced_multiplication(c[var_change_idx], a, tmp_o1_o1);
    }
    
    bitsliced_ror(tmp_o1_o1, 1);
    bitsliced_multiplication(c[var_change_idx], a, tmp_o1_o1);
    c[var_change_idx][0] ^= (c[var_change_idx][0] >> 16u);
    c[var_change_idx][1] ^= (c[var_change_idx][1] >> 16u);
    c[var_change_idx][2] ^= (c[var_change_idx][2] >> 16u);
    c[var_change_idx][3] ^= (c[var_change_idx][3] >> 16u);

    c[var_change_idx][0] &= 0xFFFF0000FFFF;
    c[var_change_idx][1] &= 0xFFFF0000FFFF;
    c[var_change_idx][2] &= 0xFFFF0000FFFF;
    c[var_change_idx][3] &= 0xFFFF0000FFFF;

    bitsliced_gf16_t tmp_o1_o1_1;
    for (unsigned int var_change_idx = (ELMTS_IN_BITSLICED / 4) - 1; var_change_idx > 0; var_change_idx--)
    {
        bitsliced_ror(tmp_o1_o1, 1);
        bitsliced_multiplication(tmp_o1_o1_1, a, tmp_o1_o1);
        bitsliced_rol_32(tmp_o1_o1_1, var_change_idx);
        bitsliced_addition(c[var_change_idx], c[var_change_idx], tmp_o1_o1_1);
    }
}
//We treat all polynomials here as if they were from the second layer

int variable_substitution_full(polynomial_t transformed, polynomial_t original, bitsliced_gf16_t variable_list[V1 + O1]){
    int i, j, k, l;
    bitsliced_gf16_t tmp, tmp2;
    memset(transformed, 0, sizeof(polynomial_t));
    int count = 0;
    int index, index1;
    // var from x_0 to x_{v_1 - 1} are transformed into x_i = x_i + \sum_{k = V1}^{V1 + O1 + O2} a_k x_k
    for (i = 0; i < V1; i++){
        index = get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES);
        bitsliced_addition(transformed[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)], 
                            transformed[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)], 
                            original[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)]);

        bitsliced_square(tmp2, variable_list[i]);
        for(j = V1; j < NUMBER_OF_VARIABLES; j++){
            index = get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES);
            index1 = get_square_index_in_polynomial(j, NUMBER_OF_VARIABLES);
            expand_bitsliced(tmp, tmp2, j - V1);
            bitsliced_muladd(transformed[get_square_index_in_polynomial(j, NUMBER_OF_VARIABLES)], 
                            original[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)], tmp);
            count++;
        }
        for(j = i + 1; j < V1; j++){
            // adding (x_i)*(x_j)
            bitsliced_addition(transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)]);

            // adding (x_i)*(a_j_l x_l) for V1 < l < N
            for (l = V1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, variable_list[j], l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(i, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
                }
            
            // adding (x_j)*(a_i_l x_l) for V1 < l < N
            for (l = V1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, variable_list[i], l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(j, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
                }
            bitsliced_gf16_t tmp_var_list[NUMB_VARS_V1_CHANGE];
            gf16_multiply_for_variable_change_v1_v1(tmp_var_list, variable_list[i], variable_list[j]);
            count+= 64;
            for(l = 0; l < NUMB_VARS_V1_CHANGE - 1; l++){
                for (  k = 0; k < ELMTS_IN_BITSLICED; k++)
                {
                    expand_bitsliced(tmp, tmp_var_list[l], k);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(k + V1, V1 + ((k + l) % ELMTS_IN_BITSLICED), 
                                    NUMBER_OF_VARIABLES)], original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], 
                                    tmp);
                    index = get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES);
                    index1 = get_binomial_index_in_polynomial(k + V1, V1 + ((k + l) % ELMTS_IN_BITSLICED), 
                            NUMBER_OF_VARIABLES);
            
                    count++;
                }
            }
            for (  k = 0; k < ELMTS_IN_BITSLICED / 2; k++)
            {
                expand_bitsliced(tmp, tmp_var_list[l], k);
                bitsliced_muladd(transformed[get_binomial_index_in_polynomial(k + V1, V1 + ((k + l) % ELMTS_IN_BITSLICED), 
                                NUMBER_OF_VARIABLES)], original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], 
                                tmp);
                index = get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES);
                index1 = get_binomial_index_in_polynomial(k + V1, V1 + ((k + l) % ELMTS_IN_BITSLICED), 
                        NUMBER_OF_VARIABLES);
        
                count++;
            }
        }
        for(j; j < V1 + O1; j++){
            // adding (x_i)*(x_j)
            bitsliced_addition(transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)]);

            // adding (x_i)*(a_j_l x_l) for V1 + O1< l < N
            for (l = V1 + O1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, variable_list[j], l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(i, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
                }
            
            // adding (x_j)*(a_i_l x_l) for V1 < l < N
            for (l = V1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, variable_list[i], l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(j, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
                }
            bitsliced_gf16_t tmp_var_list[NUMB_VARS_V1_CHANGE];
            gf16_multiply_for_variable_change_v1_o1(tmp_var_list, variable_list[i], variable_list[j]);
            count+=32;
            for(l = 0; l < ELMTS_IN_BITSLICED / 2; l++){
                for (  k = 0; k < ELMTS_IN_BITSLICED; k++)
                {
                    expand_bitsliced(tmp, tmp_var_list[l], k);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(k + V1, V1 + O1 + ((k + l) % O2), NUMBER_OF_VARIABLES)], 
                                        original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp);
                    count++;
                }
            }
        }
        for(j; j < NUMBER_OF_VARIABLES; j++){
            // adding (x_i)*(x_j)
            bitsliced_addition(transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)]);

            // adding (x_j)*(a_i_l x_l) for V1 < l < N
            for (l = V1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, variable_list[i], l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(i, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
            }
        }
    }

    for (i = V1; i < V1 + O1; i++){
        bitsliced_addition(transformed[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)], 
                            transformed[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)], 
                            original[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)]);


        bitsliced_square(tmp2, variable_list[i]);
        for(j = V1 + O1; j < NUMBER_OF_VARIABLES; j++){
            expand_bitsliced(tmp, tmp2, j - V1);
            bitsliced_muladd(transformed[get_square_index_in_polynomial(j, NUMBER_OF_VARIABLES)], 
                            original[get_square_index_in_polynomial(i, NUMBER_OF_VARIABLES)], tmp);
            count++;
        }

        for(j = i + 1; j < V1 + O1; j+=2){
            // adding (x_i)*(x_j)
            bitsliced_addition(transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)]);
            
            // If we can, we it do two by two
            if(j + 1 < V1 + O1){
                bitsliced_addition(transformed[get_binomial_index_in_polynomial(i, j + 1, NUMBER_OF_VARIABLES)],
                                transformed[get_binomial_index_in_polynomial(i, j + 1, NUMBER_OF_VARIABLES)],
                                original[get_binomial_index_in_polynomial(i, j + 1, NUMBER_OF_VARIABLES)]);
            }
            
            bitsliced_gf16_t tmp3;
            if(j + 1 < V1 + O1){
                // we can do two by two
                tmp3[0] = (variable_list[j][0] >> 32u) | (variable_list[j + 1][0] & 0xFFFFFFFF00000000u);
                tmp3[1] = (variable_list[j][1] >> 32u) | (variable_list[j + 1][1] & 0xFFFFFFFF00000000u);
                tmp3[2] = (variable_list[j][2] >> 32u) | (variable_list[j + 1][2] & 0xFFFFFFFF00000000u);
                tmp3[3] = (variable_list[j][3] >> 32u) | (variable_list[j + 1][3] & 0xFFFFFFFF00000000u);
            }else{
                tmp3[0] = (variable_list[j][0] >> 32u);
                tmp3[1] = (variable_list[j][1] >> 32u);
                tmp3[2] = (variable_list[j][2] >> 32u);
                tmp3[3] = (variable_list[j][3] >> 32u);
            }
            for (l = V1; l < NUMBER_OF_VARIABLES - O1; l++){
                expand_bitsliced(tmp2, tmp3, l - V1);
                bitsliced_muladd(transformed[get_binomial_index_in_polynomial(i, l + O1, NUMBER_OF_VARIABLES)], 
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                count++;
            }
            if(j + 1 < V1 + O1){
                for (l = V1 + O1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, tmp3, l - V1 + O1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(i, l, NUMBER_OF_VARIABLES)], 
                                original[get_binomial_index_in_polynomial(i, j + 1, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
                }
            }
        
            // adding (x_j)*(a_i_l x_l) for V1 < l < N
            for (l = V1 + O1; l < NUMBER_OF_VARIABLES; l++){
                expand_bitsliced(tmp2, tmp3, l - V1 - O1);
                bitsliced_muladd(transformed[get_binomial_index_in_polynomial(j, l, NUMBER_OF_VARIABLES)], 
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                count++;
            }
        
            // adding (x_{j+1})*(a_i_l x_l) for V1 < l < N
            if(j + 1 < V1 + O1){
                for (l = V1 + O1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, tmp3, l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(j + 1, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j + 1, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
                }
            }
            bitsliced_gf16_t c[NUMB_VARS_O1_CHANGE];
            gf16_double_multiply_for_variable_change_o1_o1(c, tmp3, variable_list[i]);
            count+=32;
            for(l = 0; l < NUMB_VARS_O1_CHANGE; l++){
                for (k = 0; k < ELMTS_IN_BITSLICED / 2; k++)
                {
                    expand_bitsliced(tmp, c[l], k);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(O1 + V1 + ((k + l) % O2), k + O1 + V1,
                                        NUMBER_OF_VARIABLES)], original[get_binomial_index_in_polynomial(i, j, 
                                        NUMBER_OF_VARIABLES)], tmp);
                    count++;
                }
                if(j + 1 < V1 + O1){
                    for (k; k < ELMTS_IN_BITSLICED; k++)
                    {
                        expand_bitsliced(tmp, c[l], k);
                        bitsliced_muladd(transformed[get_binomial_index_in_polynomial(O1 + V1 + ((k + l) % O2), k + O1 + V1,  
                                            NUMBER_OF_VARIABLES)],  original[get_binomial_index_in_polynomial(i, j + 1, 
                                            NUMBER_OF_VARIABLES)], tmp);
                        count++;
                    }
                }
            }
        }
        
        for(j = V1 + O1; j < NUMBER_OF_VARIABLES; j++){
            // adding (x_i)*(x_j)
            bitsliced_addition(transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                transformed[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)],
                                original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)]);

            // adding (x_j)*(a_i_l x_l) for V1 < l < N
            for (l = V1 + O1; l < NUMBER_OF_VARIABLES; l++){
                    expand_bitsliced(tmp2, variable_list[i], l - V1);
                    bitsliced_muladd(transformed[get_binomial_index_in_polynomial(i, l, NUMBER_OF_VARIABLES)], 
                                    original[get_binomial_index_in_polynomial(i, j, NUMBER_OF_VARIABLES)], tmp2);
                    count++;
            }
        }
    }    
    //Last O2 layer variables to put in transformed
    for(i = get_square_index_in_polynomial(O1 + V1, NUMBER_OF_VARIABLES); i < POL_NUMBER_OF_BINOMIALS; i++){
        bitsliced_addition(transformed[i], transformed[i], original[i]);
    }
    return count;
}

void polynomials_evaluation(bitsliced_gf16_t results, bitsliced_gf16_t pol[6], bitsliced_gf16_t variables_value){
    int i, j;
    bitsliced_gf16_t tmp, tmp2, tmp3;
    results[0] = 0;
    results[1] = 0;
    results[2] = 0;
    results[3] = 0;
    
    for(i = 0; i < 3; i++){
        expand_bitsliced(tmp, variables_value, i);
        bitsliced_square(tmp2, tmp);
        bitsliced_muladd(results, tmp2, pol[get_square_index_in_polynomial(i, 3)]);
        for(j = i + 1; j < 3; j++){
            expand_bitsliced(tmp2, variables_value, j);
            bitsliced_multiplication(tmp3, tmp, tmp2);
            bitsliced_muladd(results, tmp3, pol[get_binomial_index_in_polynomial(i, j, 3)]);
        }
    }
}



void polynomials_evaluation_full(bitsliced_gf16_t results, polynomial_t pol, bitsliced_gf16_t variables_value[2]){
    int i, j;
    bitsliced_gf16_t tmp, tmp2, tmp3, variables_list[NUMBER_OF_VARIABLES];
    results[0] = 0;
    results[1] = 0;
    results[2] = 0;
    results[3] = 0;
    tmp[0] = 0;
    tmp[1] = 0;
    tmp[2] = 0;
    tmp[3] = 0;
    int current_var_index = 0;
    int index_in_pol = 0;
    
    for(j = 0; j < ELMTS_IN_BITSLICED; j++){
        expand_bitsliced(variables_list[j], variables_value[0], j);
        bitsliced_muladd(tmp, variables_list[j], pol[index_in_pol]);      
        index_in_pol++;  
    }
    for(j; j < NUMBER_OF_VARIABLES; j++){
        expand_bitsliced(variables_list[j], variables_value[1], j - ELMTS_IN_BITSLICED);
        bitsliced_muladd(tmp, variables_list[j], pol[index_in_pol]);
        index_in_pol++; 
    }
    bitsliced_muladd(results, tmp, variables_list[0]);

    for(current_var_index = 1; current_var_index < NUMBER_OF_VARIABLES - 1; current_var_index++){
        bitsliced_multiplication(tmp, variables_list[current_var_index], pol[index_in_pol]);
        index_in_pol++;
        for(j = current_var_index + 1; j < NUMBER_OF_VARIABLES; j++){
            bitsliced_muladd(tmp, variables_list[j], pol[index_in_pol]);
            index_in_pol++;
        }
        bitsliced_muladd(results, tmp, variables_list[current_var_index]);
    }
    bitsliced_square(tmp, variables_list[current_var_index]);
    bitsliced_muladd(results, tmp, pol[index_in_pol]);
}

static inline uint64_t bit_parity(uint64_t v){
    v ^= v >> 32;
    v ^= v >> 16;
    v ^= v >> 8;
    v ^= v >> 4;
    v &= 0xf;
    return (0x6996 >> v) & 1;
}

void scalar_product(uint8_t pos_to_put_res, bitsliced_gf16_t res, bitsliced_gf16_t in1, bitsliced_gf16_t in2){
    uint8_t c = (pos_to_put_res & 0x3Fu);
    bitsliced_multiplication(res, in1, in2);
    res[0] = bit_parity(res[0]) << c;
    res[1] = bit_parity(res[1]) << c;
    res[2] = bit_parity(res[2]) << c;
    res[3] = bit_parity(res[3]) << c;
}

void matrix_vector_multiplication(bitsliced_gf16_t res, bitsliced_gf16_t *mat, int matrix_width, bitsliced_gf16_t vec){

    bitsliced_gf16_t tmp;
    expand_bitsliced(tmp, vec, 0);
    bitsliced_multiplication(res, tmp, mat[0]);
    for (int i = 1; i < matrix_width; i++)
    {
        expand_bitsliced(tmp, vec, i);
        bitsliced_muladd(res, tmp, mat[i]);
    }
}

void matrix_matrix_multiplication(bitsliced_gf16_t *res, 
                                bitsliced_gf16_t *mat_a, int a_width, 
                                bitsliced_gf16_t *mat_b, int b_width){

    for (int i = 0; i < b_width; i++)
    {
        matrix_vector_multiplication(res[i], mat_a, a_width, mat_b[i]);
    }    
}


// S's representation is by column
int generate_random_s(bitsliced_gf16_t s[O2]){

    uint8_t buf[4 * O2 * ((O1 + 7) / 8)];
    if(randombytes(buf, sizeof(buf)) != RNG_SUCCESS){
        return PRNG_FAILURE;
    }

    int offset = 0;
    for (size_t i = 0; i < O2; i++)
    {
        memcpy(&s[i][0], buf + offset, ((O1 + 7) / 8));
        s[i][0] &= 0xFFFFFFFF;
        offset += ((O1 + 7) / 8);
        memcpy(&s[i][1], buf + offset, ((O1 + 7) / 8));
        s[i][1] &= 0xFFFFFFFF;
        offset += ((O1 + 7) / 8);
        memcpy(&s[i][2], buf + offset, ((O1 + 7) / 8));
        s[i][2] &= 0xFFFFFFFF;
        offset += ((O1 + 7) / 8);
        memcpy(&s[i][3], buf + offset, ((O1 + 7) / 8));
        s[i][3] &= 0xFFFFFFFF;
        offset += ((O1 + 7) / 8);
    }

    return SUCCESS;

}

//T's representation is by line
int generate_random_t(bitsliced_gf16_t t[V1 + O1]){

    uint8_t buf[(V1 + O1 / 2) * sizeof(bitsliced_gf16_t)];
    if(randombytes(buf, sizeof(buf)) != RNG_SUCCESS){
        return PRNG_FAILURE;
    }

    int offset = 0;
    for (size_t i = 0; i < V1; i++)
    {
        memcpy(t[i], buf + offset, sizeof(bitsliced_gf16_t));
        offset += sizeof(bitsliced_gf16_t);
    }

    for (size_t i = V1; i < V1 + O1; i++)
    {
        memcpy(&t[i][0], buf + offset, ((O1 + 7) / 8));        
        offset += ((O1 + 7) / 8);
        memcpy(&t[i][1], buf + offset, ((O1 + 7) / 8));        
        offset += ((O1 + 7) / 8);
        memcpy(&t[i][2], buf + offset, ((O1 + 7) / 8));        
        offset += ((O1 + 7) / 8);
        memcpy(&t[i][3], buf + offset, ((O1 + 7) / 8));        
        offset += ((O1 + 7) / 8);
        shift_left_gf16(t[i], t[i], O1);
    }

    return SUCCESS;
}

//T's representation is by line
void invert_t(bitsliced_gf16_t t_inverse[V1 + O1], bitsliced_gf16_t t[V1 + O1]){

    bitsliced_gf16_t t3[O1], t1[V1], t4_transpose[V1];
    for (size_t i = 0; i < V1; i++)
    {
        t1[i][0] = t[i][0] & 0xFFFFFFFF;
        t1[i][1] = t[i][1] & 0xFFFFFFFF;
        t1[i][2] = t[i][2] & 0xFFFFFFFF;
        t1[i][3] = t[i][3] & 0xFFFFFFFF;

    }
    for (size_t i = V1; i < V1 + O1; i++)
    {
        copy_gf16(t3[i -V1], t[i]);
        copy_gf16(t_inverse[i], t[i]);
        shift_right_gf16(t3[i - V1], t3[i - V1], O1);
    }

    matrix_matrix_multiplication(t4_transpose, t3, O1, t1, V1);
    for (size_t i = 0; i < V1; i++)
    {
        shift_left_gf16(t_inverse[i], t4_transpose[i], O1);
        bitsliced_addition(t_inverse[i], t_inverse[i], t[i]);
    }    
}

//T's representation is by line
// x0, ..., x{V1-1} are in variables[0], xV1, ..., x{n-1} in variables[1]
void multiply_y_by_t(bitsliced_gf16_t result[2], bitsliced_gf16_t t[V1 + O1],
                     bitsliced_gf16_t variables[2]){

    bitsliced_gf16_t tmp, tmp1;
    copy_gf16(result[0], variables[0]);
    for (size_t i = 0; i < V1; i++)
    {
        scalar_product(i, tmp, variables[1], t[i]);
        bitsliced_addition(result[0], result[0], tmp);
    }

    copy_gf16(result[1], variables[1]);
    for (size_t i = 0; i < O1; i++)
    {
        scalar_product(i, tmp, variables[1], t[i + V1]);
        bitsliced_addition(result[1], result[1], tmp);
    }    
}

int generate_random_f(polynomial_t f){

    memset(f, 0, sizeof(polynomial_t));
    uint8_t buf[(V1*(V1 + 1)/2 + V1 * O1) * sizeof(bitsliced_gf16_t) +  (O1*(O1 + 1)/2 + O2 * O1)];
    if(randombytes(buf, sizeof(buf)) != RNG_SUCCESS){
        return PRNG_FAILURE;
    }

    int offset = 0;
    int index_in_pol = 0;
    for (size_t i = 0; i < V1; i++)
    {   
        //first and second layer
        for (size_t j = i; j < V1 + O1; j++)
        {
            memcpy(f[index_in_pol], buf + offset, sizeof(bitsliced_gf16_t));
            offset += sizeof(bitsliced_gf16_t);
            index_in_pol++;
        }
        //second layer only
        for (size_t j = V1 + O1; j < NUMBER_OF_VARIABLES; j++)
        {
            memcpy(&f[index_in_pol][0], buf + offset, O2 / 8);
            offset += O2 / 8;
            memcpy(&f[index_in_pol][1], buf + offset, O2 / 8);
            offset += O2 / 8;
            memcpy(&f[index_in_pol][2], buf + offset, O2 / 8);
            offset += O2 / 8;
            memcpy(&f[index_in_pol][3], buf + offset, O2 / 8);
            offset += O2 / 8;
            shift_left_gf16(f[index_in_pol],f[index_in_pol], O1);
            index_in_pol++;
        }        
    }
    //second layer only
    for (size_t i = V1; i < V1 + O1; i++)
    {   
        for (size_t j = i; j < NUMBER_OF_VARIABLES; j++)
        {
            memcpy(&f[index_in_pol][0], buf + offset, O1 / 8);
            offset += O1 / 8;
            memcpy(&f[index_in_pol][1], buf + offset, O1 / 8);
            offset += O1 / 8;
            memcpy(&f[index_in_pol][2], buf + offset, O1 / 8);
            offset += O1 / 8;
            memcpy(&f[index_in_pol][3], buf + offset, O1 / 8);
            offset += O1 / 8;
            shift_left_gf16(f[index_in_pol],f[index_in_pol], 32);
            index_in_pol++;
        }        
    }
    return SUCCESS;
}
