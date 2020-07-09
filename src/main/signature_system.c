#include "../lib/gf16.h"
#include "timing-functions.h"
#include <stdio.h>

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
    bitsliced_gf16_t i, j, i_times_j, result, i_square;

    i.c = 0xFFFF0000FFFF0000u;
    i.y = 0xFFFFFFFF00000000u;
    i.x = 0;
    i.y_x = 0;
    j.c = 0xAAAAAAAAAAAAAAAAu;
    j.y = 0xCCCCCCCCCCCCCCCCu;
    j.x = 0xf0f0f0f0f0f0f0f0u;
    j.y_x = 0xff00ff00ff00ff00u;

    i_times_j.c = 0x6666CCCCAAAA0000u;
    i_times_j.y = 0xAAAA6666CCCC0000u;
    i_times_j.x = 0xff0ff00f0f00000u;
    i_times_j.y_x = 0xf0f00ff0ff000000u;

    bitsliced_multiplication(&result, &i, &j);

    i.c = 0xAAAA;
    i.y = 0xCCCC;
    i.x = 0xF0F0;
    i.y_x = 0xFF00;

    i_square.c = 0x9966;
    i_square.y = 0x3C3C;
    i_square.x = 0xFF0;
    i_square.y_x = 0xFF00;

    bitsliced_square(&result, &i);


    unsigned int b;
    for (b = 0; b < 16; b++) {
        print_bitsliced(i, b);
        printf(" square = ");
        print_bitsliced(result, b);
        putchar(10);
    }

    return 0;
}
