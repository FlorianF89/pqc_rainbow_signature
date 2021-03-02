//
// Created by Florian Caullery on 7/20/20.
//

#ifndef PQC_RAINBOW_PARAMETERS_H
#define PQC_RAINBOW_PARAMETERS_H

#define O1 32
#define O2 32
#define V1 36
#define V2 (V1 + O1)
#define N (O1 + O2 + V1)
#define M (O1 + O2)
#define NUMBER_OF_EQUATIONS (O1 + O2)
#define NUMBER_OF_VARIABLES (O1 + O2 + V1)
#define POL_NUMBER_OF_BINOMIALS ((N * (N + 1)) / 2)
#define SECRET_KEY_SEED_BYTE_LENGTH 32

#endif //PQC_RAINBOW_PARAMETERS_H
