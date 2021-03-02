//
// Created by Florian Caullery on 7/20/20.
//

#include <memory.h>
#include <time.h>
#include "keygen.h"
#include "rng.h"
#include "utils_prng.h"
#include "error_codes.h"


// all expanded elements of GF(16) - GF(2)
bitsliced_gf16_t expanded_table[14] = {
    {0, -1llu, 0, 0},               //0100
    {-1llu, -1llu, 0, 0},           //1100
    {0, 0, -1llu, 0},               //0010
    {-1llu, 0, -1llu, 0},           //1010
    {0, -1llu, -1llu, 0},           //0110
    {-1llu, -1llu, -1llu, 0},       //1110
    {0, 0, 0, -1llu},               //0001
    {-1llu, 0, 0, -1llu},           //1001
    {0, -1llu, 0, -1llu},           //0101
    {-1llu, -1llu, 0, -1llu},       //1101
    {0, 0, -1llu, -1llu},           //0011
    {-1llu, 0, -1llu, -1llu},       //1011
    {0, -1llu, -1llu, -1llu},       //0111
    {-1llu, -1llu, -1llu, -1llu}    //1111
};


static inline uint64_t bit_parity(uint64_t v){
    v ^= v >> 32;
    v ^= v >> 16;
    v ^= v >> 8;
    v ^= v >> 4;
    v &= 0xf;
    return (0x6996 >> v) & 1;
}

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

void transpose_reverse_ordered_O1O2xV1_binary_matrix(uint64_t *A) {
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
    for (j = 0; j < V1; j++) {
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

static inline unsigned int get_square_index_in_polynomial(int variable_index, int number_of_variables){
    return variable_index * number_of_variables - ((variable_index * (variable_index - 1)) / 2);
}


#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
// Assumes i =< j!
static inline unsigned int get_binomial_index_in_polynomial(int i, int j, int number_of_variables){
    if(i > j) SWAP(i, j);
    return i *number_of_variables - ((i * (i - 1)) / 2) + (j - i);
}

//Assuming A[i] contains the i-th line
void print_binary_matrix(uint64_t *A, int length, int width){
    putchar(10);
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < width; j++)
        {
            printf("%ld ", (A[i] >> j) & 1u);
        }
        putchar(10);
    }
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

static void put_matrix_polynomial_into_upper_triangular_form(bitsliced_gf16_t f_times_t[NUMBER_OF_VARIABLES][2]){
    //transposing F'4 and F'7
    uint64_t A[64];
    memset(A, 0, sizeof(A));
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][0][0];
        f_times_t[i][0][0] = 0;
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < V1; i++)
    {
       f_times_t[i][1][0] ^= A[63 - i];
    }
    memset(A, 0, sizeof(A));
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][0][1];
        f_times_t[i][0][1] = 0;
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < V1; i++)
    {
       f_times_t[i][1][1] ^= A[63 - i];
    }
    memset(A, 0, sizeof(A));
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][0][2];
        f_times_t[i][0][2] = 0;
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < V1; i++)
    {
       f_times_t[i][1][2] ^= A[63 - i];
    }
    memset(A, 0, sizeof(A));
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][0][3];
        f_times_t[i][0][3] = 0;
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < V1; i++)
    {
       f_times_t[i][1][3] ^= A[63 - i];
    }

    //transposing F'5 || F'6
    //            F'8 || F'9
    // Clears unwanted values and add directly to F
    
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][1][0];
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < O1 + O2 - 1; i++)
    {

        f_times_t[i + V1][1][0] ^= (A[63 - i] & (0xFFFFFFFFFFFFFFFFllu << (i + 1)));
        f_times_t[i + V1][1][0] &= (0xFFFFFFFFFFFFFFFFllu << i);
    }
    f_times_t[O1 + O2 + V1 - 1][1][0] &= (1llu << 63);
    
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][1][1];
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < O1 + O2 - 1; i++)
    {

        f_times_t[i + V1][1][1] ^= (A[63 - i] & (0xFFFFFFFFFFFFFFFFllu << (i + 1)));
        f_times_t[i + V1][1][1] &= (0xFFFFFFFFFFFFFFFFllu << i);
    }
    f_times_t[O1 + O2 + V1 - 1][1][1] &= (1llu << 63);
    
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][1][2];
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < O1 + O2 - 1; i++)
    {

        f_times_t[i + V1][1][2] ^= (A[63 - i] & (0xFFFFFFFFFFFFFFFFllu << (i + 1)));
        f_times_t[i + V1][1][2] &= (0xFFFFFFFFFFFFFFFFllu << i);
    }
    f_times_t[O1 + O2 + V1 - 1][1][2] &= (1llu << 63);
    
    for (size_t i = V1; i < V1 + O1 + O2; i++)
    {
        A[i - V1] = f_times_t[i][1][3];
    }
    transpose_reverse_ordered_64x64_binary_matrix(A);
    for (size_t i = 0; i < O1 + O2 - 1; i++)
    {

        f_times_t[i + V1][1][3] ^= (A[63 - i] & (0xFFFFFFFFFFFFFFFFllu << (i + 1)));
        f_times_t[i + V1][1][3] &= (0xFFFFFFFFFFFFFFFFllu << i);
    }
    f_times_t[O1 + O2 + V1 - 1][1][3] &= (1llu << 63);
}


void variable_substitution_first_layer(bitsliced_gf16_t f_times_t[NUMBER_OF_VARIABLES][2], 
                                            first_layer_polynomial_t original, 
                                            bitsliced_gf16_t variable_list[V1 + O1]){
    
    memset(f_times_t, 0, 2 * NUMBER_OF_VARIABLES * sizeof(bitsliced_gf16_t));
    
    bitsliced_gf16_t tmp, tmp1;
    for (size_t i = 0; i < V1; i++)
    {
        copy_gf16(f_times_t[i][0], original[i][0]);
        copy_gf16(f_times_t[i][1], original[i][1]);
        f_times_t[i][1][0] &= 0xFFFFFFFFlu;
        f_times_t[i][1][1] &= 0xFFFFFFFFlu;
        f_times_t[i][1][2] &= 0xFFFFFFFFlu;
        f_times_t[i][1][3] &= 0xFFFFFFFFlu;
    }
    for (size_t i = 0; i < V1; i++)
    {   
        for (size_t j = 0; j < V1; j++)
        {
            expand_bitsliced(tmp, original[i][0], j);
            bitsliced_muladd(f_times_t[i][1], tmp, variable_list[j]);
        }
        for (size_t j = 0; j < O1; j++)
        {
            expand_bitsliced(tmp, original[i][1], j);
            bitsliced_muladd(f_times_t[i][1], tmp, variable_list[j + V1]);
        }
    }
    
    for (size_t i = V1; i < V1 + O1; i++)
    {   
        expand_bitsliced(tmp, variable_list[0], i - V1);
        bitsliced_multiplication(f_times_t[i][0], tmp, f_times_t[0][0]);
        bitsliced_multiplication(f_times_t[i][1], tmp, f_times_t[0][1]);
        for (size_t j = 1; j < V1; j++)
        {
            expand_bitsliced(tmp, variable_list[j], i - V1);
            bitsliced_muladd(f_times_t[i][0], tmp, f_times_t[j][0]);
            bitsliced_muladd(f_times_t[i][1], tmp, f_times_t[j][1]);
        }
    }
    for (size_t i = V1 + O1; i < V1 + O1 + O2; i++)
    {   
        expand_bitsliced(tmp, variable_list[0], i - V1);
        bitsliced_multiplication(f_times_t[i][0], tmp, f_times_t[0][0]);
        bitsliced_multiplication(f_times_t[i][1], tmp, f_times_t[0][1]);
        for (size_t j = 1; j < V1; j++)
        {
            expand_bitsliced(tmp, variable_list[j], i - V1);
            bitsliced_muladd(f_times_t[i][0], tmp, f_times_t[j][0]);
            bitsliced_muladd(f_times_t[i][1], tmp, f_times_t[j][1]);
        }
    }
    put_matrix_polynomial_into_upper_triangular_form(f_times_t);

}

void variable_substitution_second_layer(bitsliced_gf16_t f_times_t[NUMBER_OF_VARIABLES][2], 
                                            first_layer_polynomial_t original, 
                                            bitsliced_gf16_t variable_list[V1 + O1]){
    
    memset(f_times_t, 0, 2 * NUMBER_OF_VARIABLES * sizeof(bitsliced_gf16_t));
    
    bitsliced_gf16_t tmp, tmp1;
    for (size_t i = 0; i < V1 + O1; i++)
    {
        copy_gf16(f_times_t[i][0], original[i][0]);
        copy_gf16(f_times_t[i][1], original[i][1]);
    }
    // Multiply by T on the right
    for (size_t i = 0; i < V1; i++)
    {   
        for (size_t j = 0; j < V1; j++)
        {
            expand_bitsliced(tmp, original[i][0], j);
            bitsliced_muladd(f_times_t[i][1], tmp, variable_list[j]);
        }
        for (size_t j = 0; j < O1; j++)
        {
            expand_bitsliced(tmp, original[i][1], j);
            bitsliced_muladd(f_times_t[i][1], tmp, variable_list[j + V1]);
        }
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {   
        for (size_t j = V1; j < V1 + O1; j++)
        {
            expand_bitsliced(tmp, original[i][1], j- V1);
            bitsliced_muladd(f_times_t[i][1], tmp, variable_list[j]);
        }
    }    
    
    // Multiply by T^t on the left
    for (size_t i = V1 + O1; i < V1 + O1 + O2; i++)
    {   
        expand_bitsliced(tmp, variable_list[0], i - V1);
        bitsliced_multiplication(f_times_t[i][0], tmp, f_times_t[0][0]);
        bitsliced_multiplication(f_times_t[i][1], tmp, f_times_t[0][1]);
        for (size_t j = 1; j < V1; j++)
        {
            expand_bitsliced(tmp, variable_list[j], i - V1);
            bitsliced_muladd(f_times_t[i][0], tmp, f_times_t[j][0]);
            bitsliced_muladd(f_times_t[i][1], tmp, f_times_t[j][1]);
        }
        for (size_t j = V1; j < V1 + O1; j++)
        {
            expand_bitsliced(tmp, variable_list[j], i - V1);
            bitsliced_muladd(f_times_t[i][1], tmp, f_times_t[j][1]);
        }
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {   
        for (size_t j = 0; j < V1; j++)
        {
            expand_bitsliced(tmp, variable_list[j], i - V1);
            bitsliced_muladd(f_times_t[i][0], tmp, f_times_t[j][0]);
            bitsliced_muladd(f_times_t[i][1], tmp, f_times_t[j][1]);
        }
    }
    put_matrix_polynomial_into_upper_triangular_form(f_times_t);
}

//assumes the polynomial is in triangular form and constant is expanded
static void multiply_polynomial_by_constant(bitsliced_gf16_t out[NUMBER_OF_VARIABLES][2],
                                            bitsliced_gf16_t in[NUMBER_OF_VARIABLES][2],
                                            bitsliced_gf16_t constant){
    for (size_t i = 0; i < V1; i++)
    {
        bitsliced_multiplication(out[i][0], in[i][0], constant);
        bitsliced_multiplication(out[i][1], in[i][1], constant);
    }
    for (size_t i = V1; i < NUMBER_OF_VARIABLES; i++)
    {
        bitsliced_multiplication(out[i][1], in[i][1], constant);
    }
}

//adds b to a if mask is set
//assumes the polynomial is in triangular form
static void conditional_add_polynomials(bitsliced_gf16_t a[NUMBER_OF_VARIABLES][2],
                                        bitsliced_gf16_t b[NUMBER_OF_VARIABLES][2],
                                        uint64_t cond){

    for (size_t i = 0; i < V1; i++)
    {
        bitsliced_conditionnal_addition(a[i][0], a[i][0], b[i][0], cond);
        bitsliced_conditionnal_addition(a[i][1], a[i][1], b[i][1], cond);
    }
    for (size_t i = V1; i < NUMBER_OF_VARIABLES; i++)
    {
        bitsliced_conditionnal_addition(a[i][1], a[i][1], b[i][1], cond);
    }
}

int private_key_to_public_key(bitsliced_gf16_t public_key[O1 + O2][NUMBER_OF_VARIABLES][2], 
                            central_map_t original, 
                            bitsliced_gf16_t t[V1 + O1],
                            bitsliced_gf16_t s[O1]){
    int count = 0;
    for (size_t i = 0; i < O1; i++)
    {
        variable_substitution_first_layer(public_key[i], original.first_layer_polynomials[i], t);
    }
    for (size_t i = O1; i < O1 + O2; i++)
    {
        variable_substitution_second_layer(public_key[i], original.second_layer_polynomials[i - O1], t);
    }

    //Multiply by S    

    for (size_t i = O1; i < O1 + O2; i++)
    {
        uint64_t cond;
        for (size_t k = 0; k < O1; k++)
        {
            cond = bitsliced_gf16_is_one(s[k], i);
            conditional_add_polynomials(public_key[k], public_key[i], cond); 
        }
        for (size_t j = 0; j < 14; j++)
        {
            bitsliced_gf16_t f_i_times_const[NUMBER_OF_VARIABLES][2];
            memset(f_i_times_const, 0, sizeof(f_i_times_const));
            multiply_polynomial_by_constant(f_i_times_const, public_key[i], expanded_table[j]);
            bitsliced_gf16_t tmp;
            for (size_t k = 0; k < O1; k++)
            {  
                bitsliced_addition(tmp, expanded_table[j], s[k]);
                cond = gf16_is_zero(tmp, i);
                conditional_add_polynomials(public_key[k], f_i_times_const, cond); 
            }
        }
    }
    return count;
}


void polynomial_evaluation_first_layer(bitsliced_gf16_t result, first_layer_polynomial_t pol, 
                                bitsliced_gf16_t variables_value[2], int position){

    bitsliced_gf16_t accumulator[2], acc_tmp[2], tmp;

    bitsliced_multiplication(acc_tmp[0], pol[0][0], variables_value[0]);
    bitsliced_multiplication(acc_tmp[1], pol[0][1], variables_value[1]);
    expand_bitsliced(tmp, variables_value[0], 0);

    bitsliced_multiplication(accumulator[0], acc_tmp[0], tmp);
    bitsliced_multiplication(accumulator[1], acc_tmp[1], tmp);

    for (size_t i = 1; i < V1; i++)
    {
        bitsliced_multiplication(acc_tmp[0], pol[i][0], variables_value[0]);
        bitsliced_multiplication(acc_tmp[1], pol[i][1], variables_value[1]);
        expand_bitsliced(tmp, variables_value[0], i);

        bitsliced_muladd(accumulator[0], acc_tmp[0], tmp);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }

    bitsliced_addition(tmp, accumulator[0], accumulator[1]);
    result[0] = bit_parity(tmp[0]) << position;
    result[1] = bit_parity(tmp[1]) << position;
    result[2] = bit_parity(tmp[2]) << position;
    result[3] = bit_parity(tmp[3]) << position;
    
}

void polynomial_evaluation_second_layer(bitsliced_gf16_t results, first_layer_polynomial_t pol, 
                                bitsliced_gf16_t variables_value[2], int position){
    
    bitsliced_gf16_t accumulator[2], acc_tmp[2], tmp;

    bitsliced_multiplication(acc_tmp[0], pol[0][0], variables_value[0]);
    bitsliced_multiplication(acc_tmp[1], pol[0][1], variables_value[1]);
    expand_bitsliced(tmp, variables_value[0], 0);

    bitsliced_multiplication(accumulator[0], acc_tmp[0], tmp);
    bitsliced_multiplication(accumulator[1], acc_tmp[1], tmp);

    for (size_t i = 1; i < V1; i++)
    {
        bitsliced_multiplication(acc_tmp[0], pol[i][0], variables_value[0]);
        bitsliced_multiplication(acc_tmp[1], pol[i][1], variables_value[1]);
        expand_bitsliced(tmp, variables_value[0], i);

        bitsliced_muladd(accumulator[0], acc_tmp[0], tmp);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }
    for (size_t i = V1; i < V1 + O1; i++)
    {
        bitsliced_multiplication(acc_tmp[1], pol[i][1], variables_value[1]);
        expand_bitsliced(tmp, variables_value[1], i - V1);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }

    bitsliced_addition(tmp, accumulator[0], accumulator[1]);
    results[0] = bit_parity(tmp[0]) << position;
    results[1] = bit_parity(tmp[1]) << position;
    results[2] = bit_parity(tmp[2]) << position;
    results[3] = bit_parity(tmp[3]) << position;
}


void polynomial_evaluation_full(bitsliced_gf16_t results, bitsliced_gf16_t pol[NUMBER_OF_VARIABLES][2], 
                                bitsliced_gf16_t variables_value[2], int position){
    
    bitsliced_gf16_t accumulator[2], acc_tmp[2], tmp;

    bitsliced_multiplication(acc_tmp[0], pol[0][0], variables_value[0]);
    bitsliced_multiplication(acc_tmp[1], pol[0][1], variables_value[1]);
    expand_bitsliced(tmp, variables_value[0], 0);

    bitsliced_multiplication(accumulator[0], acc_tmp[0], tmp);
    bitsliced_multiplication(accumulator[1], acc_tmp[1], tmp);

    for (size_t i = 1; i < V1; i++)
    {
        bitsliced_multiplication(acc_tmp[0], pol[i][0], variables_value[0]);
        bitsliced_multiplication(acc_tmp[1], pol[i][1], variables_value[1]);
        expand_bitsliced(tmp, variables_value[0], i);

        bitsliced_muladd(accumulator[0], acc_tmp[0], tmp);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }
    for (size_t i = V1; i < NUMBER_OF_VARIABLES; i++)
    {
        bitsliced_multiplication(acc_tmp[0], pol[i][0], variables_value[0]);
        bitsliced_multiplication(acc_tmp[1], pol[i][1], variables_value[1]);
        expand_bitsliced(tmp, variables_value[1], i - V1);

        bitsliced_muladd(accumulator[0], acc_tmp[0], tmp);
        bitsliced_muladd(accumulator[1], acc_tmp[1], tmp);
    }

    bitsliced_addition(tmp, accumulator[0], accumulator[1]);
    results[0] = bit_parity(tmp[0]) << position;
    results[1] = bit_parity(tmp[1]) << position;
    results[2] = bit_parity(tmp[2]) << position;
    results[3] = bit_parity(tmp[3]) << position;
    
}

void scalar_product(uint8_t pos_to_put_res, bitsliced_gf16_t res, bitsliced_gf16_t in1, bitsliced_gf16_t in2){
    uint8_t c = (pos_to_put_res & 0x3Fu);
    bitsliced_multiplication(res, in1, in2);
    res[0] = bit_parity(res[0]) << c;
    res[1] = bit_parity(res[1]) << c;
    res[2] = bit_parity(res[2]) << c;
    res[3] = bit_parity(res[3]) << c;
}

void column_matrix_vector_multiplication(bitsliced_gf16_t res, bitsliced_gf16_t *mat, int matrix_width, bitsliced_gf16_t vec){

    bitsliced_gf16_t tmp;
    expand_bitsliced(tmp, vec, 0);
    bitsliced_multiplication(res, tmp, mat[0]);
    for (int i = 1; i < matrix_width; i++)
    {
        expand_bitsliced(tmp, vec, i);
        bitsliced_muladd(res, tmp, mat[i]);
    }
}

void column_matrices_multiplication(bitsliced_gf16_t *res, 
                                bitsliced_gf16_t *mat_a, int a_width, 
                                bitsliced_gf16_t *mat_b, int b_width){

    for (int i = 0; i < b_width; i++)
    {
        column_matrix_vector_multiplication(res[i], mat_a, a_width, mat_b[i]);
    }    
}

void line_matrix_vector_multiplication(bitsliced_gf16_t res, 
                                bitsliced_gf16_t *mat, 
                                int matrix_length, 
                                bitsliced_gf16_t vec){

    scalar_product(0, res, mat[0], vec);
    for (int i = 1; i < matrix_length; i++)
    {
        bitsliced_gf16_t tmp;
        scalar_product(i, tmp, mat[i], vec);
        bitsliced_addition(res, tmp, res);
    }
    
}

void line_matrices_multiplication(bitsliced_gf16_t *res, 
                                bitsliced_gf16_t *mat_a, int a_length,
                                bitsliced_gf16_t *mat_b, int b_length){
    
    bitsliced_gf16_t tmp;
    for (int i = 0; i < a_length; i++)
    {
        expand_bitsliced(tmp, mat_a[i], 0);
        bitsliced_multiplication(res[i], tmp, mat_b[0]);
        for (int j = 1; j < b_length; j++)
        {
            expand_bitsliced(tmp, mat_a[i], j);
            bitsliced_muladd(res[i], tmp, mat_b[j]);
        }
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
        s[i][0] <<= 32;
        offset += ((O1 + 7) / 8);
        memcpy(&s[i][1], buf + offset, ((O1 + 7) / 8));
        s[i][1] <<= 32;
        offset += ((O1 + 7) / 8);
        memcpy(&s[i][2], buf + offset, ((O1 + 7) / 8));
        s[i][2] <<= 32;
        offset += ((O1 + 7) / 8);
        memcpy(&s[i][3], buf + offset, ((O1 + 7) / 8));
        s[i][3] <<= 32;
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

    column_matrices_multiplication(t4_transpose, t3, O1, t1, V1);
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

//S's representation is by line, only the second half of the bitsliced elements is written in S
void multiply_y_by_s(bitsliced_gf16_t result, bitsliced_gf16_t s[O1],
                     bitsliced_gf16_t variables){

    bitsliced_gf16_t tmp, tmp1;
    copy_gf16(result, variables);
    for (size_t i = 0; i < O1; i++)
    {
        scalar_product(i, tmp, variables, s[i]);
        bitsliced_addition(result, result, tmp);
    }   
}

int generate_random_central_map(central_map_t *f){

    for (size_t j = 0; j < O1; j++)
    {    
        memset(f->first_layer_polynomials[j], 0, sizeof(first_layer_polynomial_t));
        for (size_t i = 0; i < V1; i++)
        {
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][0][0], (V1 + 7) / 8);
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][0][1], (V1 + 7) / 8);
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][0][2], (V1 + 7) / 8);
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][0][3], (V1 + 7) / 8);
            f->first_layer_polynomials[j][i][0][0] = (f->first_layer_polynomials[j][i][0][0] << i) & (1llu << (V1)) - 1;
            f->first_layer_polynomials[j][i][0][1] = (f->first_layer_polynomials[j][i][0][1] << i) & (1llu << (V1)) - 1;
            f->first_layer_polynomials[j][i][0][2] = (f->first_layer_polynomials[j][i][0][2] << i) & (1llu << (V1)) - 1;
            f->first_layer_polynomials[j][i][0][3] = (f->first_layer_polynomials[j][i][0][3] << i) & (1llu << (V1)) - 1;
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][1][0], (O1 + 7) / 8);
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][1][1], (O1 + 7) / 8);
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][1][2], (O1 + 7) / 8);
            randombytes((uint8_t *) &f->first_layer_polynomials[j][i][1][3], (O1 + 7) / 8);
            f->first_layer_polynomials[j][i][1][0] &= (1llu << O1) - 1;
            f->first_layer_polynomials[j][i][1][1] &= (1llu << O1) - 1;
            f->first_layer_polynomials[j][i][1][2] &= (1llu << O1) - 1;
            f->first_layer_polynomials[j][i][1][3] &= (1llu << O1) - 1;
        }
    }
    for (size_t j = 0; j < O2; j++)
    {
        memset(f->second_layer_polynomials[j], 0, sizeof(second_layer_polynomial_t));
        for (size_t i = 0; i < V1; i++)
        {
            randombytes((uint8_t *) &f->second_layer_polynomials[j][i][0][0], (V1 + 7) / 8);
            randombytes((uint8_t *) &f->second_layer_polynomials[j][i][0][1], (V1 + 7) / 8);
            randombytes((uint8_t *) &f->second_layer_polynomials[j][i][0][2], (V1 + 7) / 8);
            randombytes((uint8_t *) &f->second_layer_polynomials[j][i][0][3], (V1 + 7) / 8);
            randombytes((uint8_t *) f->second_layer_polynomials[j][i][1], sizeof(bitsliced_gf16_t));
        }
        for (size_t i = V1; i < V1 + O1; i++)
        {
            randombytes((uint8_t *) f->second_layer_polynomials[j][i][1], sizeof(bitsliced_gf16_t));        
        }
    }
    return SUCCESS;
}