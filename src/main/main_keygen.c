#include <stdio.h>
#include <string.h>
#include "../lib/gf16.h"
#include "../lib/parameters.h"
#include "../lib/keygen.h"
#include "../lib/sign.h"
#include "../lib/verify.h"
#include "../lib/rng.h"
#include "../lib/error_codes.h"
#include "../lib/utils_hash.h"
#include <stdlib.h>


int byte_fdump(FILE * fp, const unsigned char *v, unsigned n_byte)
{
	int r = 0;

	for(unsigned i=0;i<n_byte;i++) {
		int t = fprintf( fp , "%02x", v[i] );
		if( 0 > t ) return t;
		else r += t;
	}
	return r;
}

int main( int argc , char ** argv )
{
    printf( "Rainbow Classic Ia\n");
    printf("sk size: %lu\n", 48 );
    printf("pk size: %lu\n",  NUMBER_OF_EQUATIONS * NUMBER_OF_VARIABLES * 2 * sizeof(bitsliced_gf16_t));
    printf("hash size: %d\n", _HASH_LEN );
    printf("signature size: %lu\n\n", 2 * sizeof(bitsliced_gf16_t));

	if( !((3 == argc) || (4 == argc)) ) {
		printf("Usage:\n\n\trainbow-genkey pk_file_name sk_file_name [random_seed_file]\n\n");
		return -1;
	}

    uint8_t seed[48] = {0};
    randombytes_init(seed, NULL, NULL);
    
    central_map_t f;
    bitsliced_gf16_t public_key[O1 + O2][NUMBER_OF_VARIABLES][2], s[O1], t[O1 + V1], 
                tmp, tmp1, tmp2, variables[2], var_times_t[2];;

    generate_random_central_map(&f);    
    generate_random_t(t);    
    generate_random_s(s);

    private_key_to_public_key(public_key, f, t, s);

    return 0;
}