/// @file rainbow_keypair.c
/// @brief implementations of functions in rainbow_keypair.h
///
#include "rainbow_keypair.h"
#include "rainbow_keypair_computation.h"

#include "blas_comm.h"
#include "blas.h"
#include "rainbow_blas.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>


/////////////////////////////////////////////////////////////////


#include "utils_prng.h"


static
void generate_S_T( unsigned char * s_and_t , prng_t * prng0 )
{
    prng_gen( prng0 , s_and_t , _O1_BYTE*_O2 ); // S1
    s_and_t += _O1_BYTE*_O2;
    prng_gen( prng0 , s_and_t , _V1_BYTE*_O1 ); // T1
    s_and_t += _V1_BYTE*_O1;
    prng_gen( prng0 , s_and_t , _V1_BYTE*_O2 ); // T2
    s_and_t += _V1_BYTE*_O2;
    prng_gen( prng0 , s_and_t , _O1_BYTE*_O2 ); // T3
}


static
unsigned generate_l1_F12( unsigned char * sk, prng_t * prng0 )
{
    unsigned n_byte_generated = 0;
    prng_gen( prng0 , sk , _O1_BYTE * N_TRIANGLE_TERMS(_V1) ); // l1_F1
    sk += _O1_BYTE * N_TRIANGLE_TERMS(_V1);
    n_byte_generated += _O1_BYTE * N_TRIANGLE_TERMS(_V1);

    prng_gen( prng0 , sk , _O1_BYTE * _V1*_O1 );  // l1_F2
    sk += _O1_BYTE * _V1*_O1;
    n_byte_generated += _O1_BYTE * _V1*_O1;

    return n_byte_generated;
}


static
unsigned generate_l2_F12356( unsigned char * sk, prng_t * prng0 )
{
    unsigned n_byte_generated = 0;

    prng_gen( prng0 , sk , _O2_BYTE * N_TRIANGLE_TERMS(_V1) ); // l2_F1
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_V1);
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_V1);

    prng_gen( prng0 , sk , _O2_BYTE * _V1*_O1 ); // l2_F2
    sk += _O2_BYTE * _V1*_O1;
    n_byte_generated += _O2_BYTE * _V1*_O1;

    prng_gen( prng0 , sk , _O2_BYTE * _V1*_O2 ); // l2_F3
    sk += _O2_BYTE * _V1*_O1;
    n_byte_generated += _O2_BYTE * _V1*_O1;

    prng_gen( prng0 , sk , _O2_BYTE * N_TRIANGLE_TERMS(_O1) ); // l2_F5
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_O1);
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_O1);

    prng_gen( prng0 , sk , _O2_BYTE * _O1*_O2 ); // l2_F6
    n_byte_generated += _O2_BYTE * _O1*_O2;

    return n_byte_generated;
}


static
void generate_B1_B2( unsigned char * sk , prng_t * prng0 )
{
    sk += generate_l1_F12( sk , prng0 );
    generate_l2_F12356( sk , prng0 );
}


//////////////////////////////////////////////////////////



void cpk_to_pk( pk_t * rpk, const cpk_t * cpk )
{
    // procedure:  cpk_t --> extcpk_t  --> pk_t

    // convert from cpk_t to extcpk_t
    ext_cpk_t * pk = (ext_cpk_t *) aligned_alloc( 32, sizeof(ext_cpk_t) );
    // setup prng
    prng_t prng0;
    prng_set( &prng0 , cpk->pk_seed , LEN_SKSEED );

    // generating parts of key with prng
    generate_l1_F12( pk->l1_Q1 , &prng0 );
    // copying parts of key from input. l1_Q3, l1_Q5, l1_Q6, l1_Q9
    memcpy( pk->l1_Q3 , cpk->l1_Q3 , _O1_BYTE*( _V1*_O2 + N_TRIANGLE_TERMS(_O1) + _O1*_O2 + N_TRIANGLE_TERMS(_O2) ) );

    // generating parts of key with prng
    generate_l2_F12356( pk->l2_Q1 , &prng0 );
    // copying parts of key from input: l2_Q9
    memcpy( pk->l2_Q9 , cpk->l2_Q9 , _O2_BYTE* N_TRIANGLE_TERMS(_O2) );

    // convert from extcpk_t to pk_t
    extcpk_to_pk( rpk , pk );

    free( pk );
}



/////////////////////////////////////////////////////////



static
void calculate_t4( unsigned char * t2_to_t4 , const unsigned char *t1 , const unsigned char *t3 )
{
    //  t4 = T_sk.t1 * T_sk.t3 - T_sk.t2
    unsigned char temp[_V1_BYTE+32];
    unsigned char * t4 = t2_to_t4;
    for(unsigned i=0;i<_O2;i++) {  /// t3 width
        gfmat_prod( temp , t1 , _V1_BYTE , _O1 , t3 );
        gf256v_add( t4 , temp , _V1_BYTE );
        t4 += _V1_BYTE;
        t3 += _O1_BYTE;
    }
}



static
void obsfucate_l1_polys( unsigned char * l1_polys , const unsigned char * l2_polys , unsigned n_terms , const unsigned char * s1 )
{
    unsigned char temp[_O1_BYTE + 32];
    while( n_terms-- ) {
        gfmat_prod( temp , s1 , _O1_BYTE , _O2 , l2_polys );
        gf256v_add( l1_polys , temp , _O1_BYTE );
        l1_polys += _O1_BYTE;
        l2_polys += _O2_BYTE;
    }
}



///////////////////  Classic //////////////////////////////////


static
void _generate_secretkey( sk_t* sk, const unsigned char *sk_seed )
{
    memcpy( sk->sk_seed , sk_seed , LEN_SKSEED );

    // set up prng
    prng_t prng0;
    prng_set( &prng0 , sk_seed , LEN_SKSEED );

    // generating secret key with prng.
    generate_S_T( sk->s1 , &prng0 );
    generate_B1_B2( sk->l1_F1 , &prng0 );

    // clean prng
    memset( &prng0 , 0 , sizeof(prng_t) );
}


void generate_secretkey( sk_t* sk, const unsigned char *sk_seed )
{
    _generate_secretkey( sk , sk_seed );
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );
}



void generate_keypair( pk_t * rpk, sk_t* sk, const unsigned char *sk_seed )
{
    _generate_secretkey( sk , sk_seed );

    // set up a temporary structure ext_cpk_t for calculating public key.
    ext_cpk_t * pk = (ext_cpk_t*) aligned_alloc( 32 , sizeof(ext_cpk_t) );
    calculate_Q_from_F( pk, sk , sk );   // compute the public key in ext_cpk_t format.
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );

    obsfucate_l1_polys( pk->l1_Q1 , pk->l2_Q1 , N_TRIANGLE_TERMS(_V1) , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q2 , pk->l2_Q2 , _V1*_O1 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q3 , pk->l2_Q3 , _V1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q5 , pk->l2_Q5 , N_TRIANGLE_TERMS(_O1) , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q6 , pk->l2_Q6 , _O1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q9 , pk->l2_Q9 , N_TRIANGLE_TERMS(_O2) , sk->s1 );
    // so far, the pk contains the full pk but in ext_cpk_t format.

    extcpk_to_pk( rpk , pk );     // convert the public key from ext_cpk_t to pk_t.
    free( pk );
}



/////////////////////   Cyclic   //////////////////////////////////


void generate_secretkey_cyclic( sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed )
{
    memcpy( sk->sk_seed , sk_seed , LEN_SKSEED );

    // prng for sk
    prng_t prng0;
    prng_set( &prng0 , sk_seed , LEN_SKSEED );
    generate_S_T( sk->s1 , &prng0 );
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );

    // prng for pk
    sk_t inst_Qs;
    sk_t * Qs = &inst_Qs;
    prng_t prng1;
    prng_set( &prng1 , pk_seed , LEN_PKSEED );
    generate_B1_B2( Qs->l1_F1 , &prng1 );

    obsfucate_l1_polys( Qs->l1_F1 , Qs->l2_F1 , N_TRIANGLE_TERMS(_V1) , sk->s1 );
    obsfucate_l1_polys( Qs->l1_F2 , Qs->l2_F2 , _V1*_O1 , sk->s1 );

    // calcuate the parts of sk according to pk.
    calculate_F_from_Q( sk , Qs , sk );

    // clean prng for sk
    memset( &prng0 , 0 , sizeof(prng_t) );
}





void generate_keypair_cyclic( cpk_t * pk, sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed )
{
    memcpy( pk->pk_seed , pk_seed , LEN_PKSEED );

    // prng for sk
    prng_t prng;
    prng_t * prng0 = &prng;
    prng_set( prng0 , sk_seed , LEN_SKSEED );
    generate_S_T( sk->s1 , prng0 );   // S,T:  only a part of sk

    unsigned char * t2 = (unsigned char *) aligned_alloc( 32, sizeof(sk->t4) );
    memcpy( t2 , sk->t4 , _V1_BYTE*_O2 );        // temporarily store t2
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );    // t2 <- t4

    // prng for pk
    sk_t inst_Qs;
    sk_t * Qs = &inst_Qs;
    prng_t * prng1 = &prng;
    prng_set( prng1 , pk_seed , LEN_PKSEED );
    generate_B1_B2( Qs->l1_F1 , prng1 );  // generating l1_Q1, l1_Q2, l2_Q1, l2_Q2, l2_Q3, l2_Q5, l2_Q6
    obsfucate_l1_polys( Qs->l1_F1 , Qs->l2_F1 , N_TRIANGLE_TERMS(_V1) , sk->s1 );
    obsfucate_l1_polys( Qs->l1_F2 , Qs->l2_F2 , _V1*_O1 , sk->s1 );
    // so far, the Qs contains l1_F1, l1_F2, l2_F1, l2_F2, l2_F3, l2_F5, l2_F6.

    calculate_F_from_Q( sk , Qs , sk );          // calcuate the rest parts of secret key from Qs and S,T

    memcpy( sk->t4 , t2 , _V1_BYTE*_O2 );        // restore t2
    calculate_Q_from_F_cyclic( pk, sk , sk );    // calculate the rest parts of public key: l1_Q3, l1_Q5, l1_Q6, l1_Q9, l2_Q9

    obsfucate_l1_polys( pk->l1_Q3 , Qs->l2_F3 , _V1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q5 , Qs->l2_F5 , N_TRIANGLE_TERMS(_O1) , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q6 , Qs->l2_F6 , _O1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q9 , pk->l2_Q9 , N_TRIANGLE_TERMS(_O2) , sk->s1 );

    // clean
    memset( &prng , 0 , sizeof(prng_t) );
    memset( t2 , 0 , sizeof(sk->t4) );
    free( t2 );
}



void generate_compact_keypair_cyclic( cpk_t * pk, csk_t* rsk, const unsigned char *pk_seed , const unsigned char *sk_seed )
{
    memcpy( rsk->pk_seed , pk_seed , LEN_PKSEED );
    memcpy( rsk->sk_seed , sk_seed , LEN_SKSEED );

    sk_t * sk = (sk_t *) aligned_alloc( 32 , sizeof(sk_t) );
    generate_keypair_cyclic( pk , sk , pk_seed , sk_seed );
    memset( sk , 0 , sizeof(sk_t) );
    free( sk );    // dispose of sk. don't need to output.
}




