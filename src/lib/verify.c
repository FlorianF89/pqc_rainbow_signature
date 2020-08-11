//
// Created by Florian Caullery on 8/10/20.
//

#include <string.h>
#include "verify.h"


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
        copy_gf16(&tmp, &x0_x95_expanded[i]);
        bitsliced_square(&tmp1, &tmp);
        bitsliced_multiplication(&tmp2, &tmp1, &f->coefficients[offset]);
        bitsliced_addition(evaluations, evaluations, &tmp2);
        offset++;
        for (j = 1; j < N - i; j++) {
            bitsliced_multiplication(&tmp1, &tmp, &f->coefficients[offset]);
            bitsliced_multiplication(&tmp2, &tmp1, &x0_x95_expanded[j]);
            bitsliced_addition(evaluations, evaluations, &tmp2);
            offset++;
        }
    }
}