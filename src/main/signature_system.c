#include "../lib/gf16.h"
#include "timing-functions.h"
#include "../lib/keygen.h"
#include "../lib/sign.h"
#include "../lib/verify.h"
#include <stdio.h>
#include <x86intrin.h>
#include <string.h>
#include <time.h>

static void print_bitsliced(const bitsliced_gf16_t in, unsigned int bit_position) {
    int is_first = 0;
    if ((in.c >> bit_position) & 0x01u) {
        printf("1");
        is_first = 1;
    }
    if ((in.y >> bit_position) & 0x01u) {
        if (is_first) {
            printf(" + ");
        }
        printf("y");
        is_first = 1;
    }
    if ((in.x >> bit_position) & 0x01u) {
        if (is_first) {
            printf(" + ");
        }
        printf("X");
        is_first = 1;
    }
    if ((in.y_x >> bit_position) & 0x01u) {
        if (is_first) {
            printf(" + ");
        }
        printf("y*X");
        is_first = 1;
    }
    if (is_first == 0) {
        printf("0");
    }
}

int main(int argc, char **argv) {
    private_key_t private_key;
    memset(&private_key, 0x00, sizeof(private_key));
    public_key_t public_key;
    memset(&public_key, 0x00, sizeof(public_key));
    uint8_t seed[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_t prng;
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    int i, j;
    int total_iterations = 100;
    unsigned int dummy;
    unsigned long long t1, t2, total = 0;
    clock_t begin = clock();
    for (j = 0; j < total_iterations; j++) {
        generate_private_key(&private_key, &prng);
//        t1 = __rdtscp(&dummy);
        derive_public_key_from_private_key(&public_key, &private_key);
//        t2 = __rdtscp(&dummy);
//        total += t2 - t1;
    }
    clock_t end = clock();
    clock_t time_spent = end - begin;
    printf("generate_keypair = %.2f\n", (double) time_spent / total_iterations);

    bitsliced_gf16_t signature[2];
    memset(signature, 0x00, sizeof(signature));

    uint8_t message[SECRET_KEY_SEED_BYTE_LENGTH];
    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    memset(message, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);


    total = 0;
    begin = clock();
    for (i = 0; i < total_iterations; i++) {
//        t1 = __rdtscp(&dummy);
        rainbow_sign(signature, &private_key, message, &prng, SECRET_KEY_SEED_BYTE_LENGTH);
//        t2 = __rdtscp(&dummy);
//        total += t2 - t1;
    }
    end = clock();
    time_spent = end - begin;
    printf("generate_keypair = %.2f\n", (double) time_spent / total_iterations);
//    printf("full sign: %.2fcc\n", (double) total / total_iterations);

    bitsliced_gf16_t x0_x63, x64_x95, evaluation;
    bitsliced_quadratic_polynomials_t f;

    memset(seed, 0x00, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_set(&prng, seed, SECRET_KEY_SEED_BYTE_LENGTH);
    prng_gen(&prng, (uint8_t *) &x0_x63, sizeof(bitsliced_gf16_t));

    for (i = 0; i < N * (N + 1) / 2; i++) {
        prng_gen(&prng, (uint8_t *) &f.coefficients[i], sizeof(bitsliced_gf16_t));
    }

    total = 0;
    begin = clock();
    for (i = 0; i < total_iterations; i++) {
//        t1 = __rdtscp(&dummy);
        evaluate_quadratic_polynomials(&evaluation, &f, &x0_x63, &x64_x95);
//        t2 = __rdtscp(&dummy);
//        total += t2 - t1;
    }
    end = clock();
    time_spent = end - begin;
    printf("generate_keypair = %.2f\n", (double) time_spent / total_iterations);
//    printf("verify: %.2fcc\n", (double) (total) / total_iterations);
    return 0;
}
