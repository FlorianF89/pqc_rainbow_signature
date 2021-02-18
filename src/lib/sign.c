//
// Created by Florian Caullery on 8/2/20.
//

#include <memory.h>
#include "sign.h"
#include "keygen.h"
#include "hash_len_config.h"
#include "utils_hash.h"

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


void print_32x32_gf16_system(bitsliced_gf16_t m[32], bitsliced_gf16_t solution) {
    int i, j;
    for (i = 0; i < 32; i++) {
        printf("l%02d: ", i);
        for (j = 0; j < 32; j++) {
            print_bitsliced(m[j], i);
        }
        print_bitsliced(solution, i);
        putchar(10);
    }
    putchar(10);
}

 void swap_elements_i_and_k(bitsliced_gf16_t r, bitsliced_gf16_t in, unsigned int i, unsigned int j) {

    uint64_t b = in[0];
    uint64_t x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r[0] = b ^ ((x << i) | (x << j));

    b = in[1];
    x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r[1] = b ^ ((x << i) | (x << j));

    b = in[2];
    x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r[2] = b ^ ((x << i) | (x << j));

    b = in[3];
    x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
    r[3] = b ^ ((x << i) | (x << j));

}


static inline void expand_bitsliced(bitsliced_gf16_t out, bitsliced_gf16_t const in, const unsigned int pos){
    
    uint8_t c = (pos & 0x3Fu);
    out[0] = 0 - ((in[0] >> c) & 0x1llu);
    out[1] = 0 - ((in[1] >> c) & 0x1llu);
    out[2] = 0 - ((in[2] >> c) & 0x1llu);
    out[3] = 0 - ((in[3] >> c) & 0x1llu);
}

//todo Does not work if system is undetermined
int solve_32x32_gf16_system(bitsliced_gf16_t solution, bitsliced_gf16_t equations_coefficients[32],
                            bitsliced_gf16_t linear_coefficients) {
    copy_gf16(solution, linear_coefficients);
    solution[0] &= 0xFFFFFFFF;
    solution[1] &= 0xFFFFFFFF;
    solution[2] &= 0xFFFFFFFF;
    solution[3] &= 0xFFFFFFFF;
    unsigned int i, j, k;
    bitsliced_gf16_t tmp, tmp1, tmp2;
    for (i = 0; i < 31; i++) {
        for (j = i + 1; j < 31; j++) {
            if (gf16_is_zero(equations_coefficients[i], i)) {
//                printf("\n swapping lines %d and %d. Before:\n", i, j);
//                print_32x32_gf16_system(equations_coefficients, *solution);
                for (k = i; k < 32; k++)
                    swap_elements_i_and_k(equations_coefficients[k], equations_coefficients[k], i, j);
                swap_elements_i_and_k(solution, solution, i, j);
//                printf("\n swapping lines %d and %d. After:\n", i, j);
//                print_32x32_gf16_system(equations_coefficients, *solution);
            } else {
                break;
            }
        }
        if (j == 32) {
            return 0;
        }
        expand_bitsliced(tmp, equations_coefficients[i], i);
        shift_left_gf16(tmp, tmp, i + 1);
        bitsliced_inversion(tmp1, tmp);
        bitsliced_multiplication(tmp, tmp1, equations_coefficients[i]);
        equations_coefficients[i][0] &= ((((uint_fast64_t) 1) << (i + 1)) - 1);
        equations_coefficients[i][1] &= ((((uint_fast64_t) 1) << (i + 1)) - 1);
        equations_coefficients[i][2] &= ((((uint_fast64_t) 1) << (i + 1)) - 1);
        equations_coefficients[i][3] &= ((((uint_fast64_t) 1) << (i + 1)) - 1);
        for (j = i + 1; j < 32; j++)
        {
            expand_bitsliced(tmp1, equations_coefficients[j], i);
            bitsliced_muladd(equations_coefficients[j], tmp1, tmp);
        }
        expand_bitsliced(tmp1, solution, i);
        bitsliced_muladd(solution, tmp1, tmp);
    }

    for (i = 31; ((int) i) > 0; i--){
        
        expand_bitsliced(tmp, equations_coefficients[i], i);
        bitsliced_inversion(tmp1, tmp); //tmp1 = (1/m_{i,i},)
        shift_right_gf16(tmp1, tmp1, 31 - i);

        expand_bitsliced(tmp, solution, i);
        shift_right_gf16(tmp, tmp, 31 - i);
        bitsliced_multiplication(tmp2, tmp, tmp1); 

        tmp[0] = (uint32_t) equations_coefficients[i][0] & (( 1 << i) - 1) | (1 << i);
        tmp[1] = (uint32_t) equations_coefficients[i][1] & (( 1 << i) - 1);
        tmp[2] = (uint32_t) equations_coefficients[i][2] & (( 1 << i) - 1);
        tmp[3] = (uint32_t) equations_coefficients[i][3] & (( 1 << i) - 1);

        solution[0] &= (uint32_t) (~(1 << i));
        solution[1] &= (uint32_t) (~(1 << i));
        solution[2] &= (uint32_t) (~(1 << i));
        solution[3] &= (uint32_t) (~(1 << i));

        bitsliced_muladd(solution, tmp, tmp2);
        equations_coefficients[i][0] = 1 << i;
        equations_coefficients[i][1] = 0;
        equations_coefficients[i][2] = 0;
        equations_coefficients[i][3] = 0;
    }
    bitsliced_inversion(tmp, equations_coefficients[0]);
    tmp[0] |= 0xFFFFFFFE;
    bitsliced_multiplication(tmp1, tmp, solution);
    copy_gf16(solution, tmp1);
    return 1;
}

void polynomials_evaluation_in_v1_variables(bitsliced_gf16_t results, polynomial_t pol, bitsliced_gf16_t variables_value){
    int i, j;
    bitsliced_gf16_t tmp, tmp2, tmp3, variables_list[V1];
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
    
    for(j = 0; j < V1; j++){
        expand_bitsliced(variables_list[j], variables_value, j);
        bitsliced_muladd(tmp, variables_list[j], pol[index_in_pol]);      
        index_in_pol++;  
    }
    index_in_pol += NUMBER_OF_VARIABLES - V1;
    bitsliced_muladd(results, tmp, variables_list[0]);

    for(current_var_index = 1; current_var_index < V1 - 1; current_var_index++){
        bitsliced_multiplication(tmp, variables_list[current_var_index], pol[index_in_pol]);
        index_in_pol++;
        for(j = current_var_index + 1; current_var_index < V1; current_var_index++){
            bitsliced_muladd(tmp, variables_list[j], pol[index_in_pol]);
            index_in_pol++;
        }
        bitsliced_muladd(results, tmp, variables_list[current_var_index]);
        index_in_pol += NUMBER_OF_VARIABLES - V1;
    }
    bitsliced_square(tmp, variables_list[current_var_index]);
    bitsliced_muladd(results, tmp, pol[index_in_pol]);
}

void prepare_first_layer_system(bitsliced_gf16_t system[O1], polynomial_t pol, bitsliced_gf16_t variable_list[V1]){

    for (size_t i = 0; i < O1; i++)
    {
        int index_in_pol = V1 + i;
        bitsliced_multiplication(system[i], pol[index_in_pol], variable_list[0]);
        for (size_t j = 1; j < V1; j++)
        {
            index_in_pol += NUMBER_OF_VARIABLES  - j;
            bitsliced_muladd(system[i], variable_list[j], pol[index_in_pol]);
        }
    }
}

void polynomials_evaluation_in_o1_variables(bitsliced_gf16_t results, polynomial_t pol, 
                                            bitsliced_gf16_t variables_v1_values,
                                            bitsliced_gf16_t variables_o1_values,
                                            bitsliced_gf16_t previous_evaluation_in_v1){
    int i, j;
    bitsliced_gf16_t tmp, tmp2, tmp3, variables_list[O1];
    results[0] = previous_evaluation_in_v1[0];
    results[1] = previous_evaluation_in_v1[1];
    results[2] = previous_evaluation_in_v1[2];
    results[3] = previous_evaluation_in_v1[3];
    int current_var_index = 0;
    int index_in_pol = V1;
    
    expand_bitsliced(variables_list[0], variables_o1_values, 0);
    bitsliced_multiplication(tmp, pol[index_in_pol], variables_list[0]);
    for (size_t i = 1; i < O1; i++)
    {   
        index_in_pol++;
        expand_bitsliced(variables_list[i], variables_o1_values, i);
        bitsliced_muladd(tmp, variables_list[i], pol[index_in_pol]);
    }

    expand_bitsliced(tmp2, variables_v1_values, 0);
    bitsliced_muladd(results, tmp2, tmp);

    for (size_t i = 1; i < V1; i++)
    {
        index_in_pol += O2 + V1 - i;
        bitsliced_multiplication(tmp, pol[index_in_pol], variables_list[0]);
        for (size_t j = 0; j < O1; j++)
        {
            index_in_pol++;
            bitsliced_muladd(tmp, variables_list[j], pol[index_in_pol]);
        }
        expand_bitsliced(tmp2, variables_v1_values, i);
        bitsliced_muladd(results, tmp2, tmp);
    }
}



void prepare_second_layer_system(bitsliced_gf16_t system[O2], polynomial_t pol, bitsliced_gf16_t variable_list[O1 + V1]){

    for (size_t i = 0; i < O2; i++){
        int index_in_pol = V1 + O1 + i;
        bitsliced_multiplication(system[i], pol[index_in_pol], variable_list[0]);
        for (size_t j = 1; j < V1 + O1; j++){
            index_in_pol += NUMBER_OF_VARIABLES  - j;
            bitsliced_muladd(system[i], variable_list[j], pol[index_in_pol]);
        }
    }
}
