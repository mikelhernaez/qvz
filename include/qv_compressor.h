//
//  fasta_compressor.h
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef iDoComp_v1_fasta_compressor_h

#define iDoComp_v1_fasta_compressor_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CARDINALITY 1000000
#define MAX_ALPHA 300000000

#define m_INTS  18
#define m_CHARS 13
#define m_SIGNS 23

//#define mI 18
//#define mc 13
//#define mS 23

#define COMPRESSION 0
#define DECOMPRESSION 1

typedef enum BASEPAIRS{bp_A,bp_C,bp_G,bp_T,bp_N,bp_U,bp_R,bp_Y,bp_K,bp_M,bp_S,bp_W,bp_B,bp_D,bp_H, bp_V, bp_X, bp_a,bp_c,bp_g,bp_t,bp_n,bp_u,bp_r,bp_y,bp_k,bp_m,bp_s,bp_w,bp_b,bp_d,bp_h, bp_v, bp_x} BASEPAIR;

typedef struct Arithmetic_code_t{
    
    uint32_t scale3;
    uint32_t m;
    uint32_t l;
    uint32_t u;
    uint32_t t;
}*Arithmetic_code;

typedef struct osStream_t{
    
    uint8_t bufferByte;
    int8_t bitPos;
    FILE *fos;
    uint32_t osLength;
    uint8_t *osBuffer;
    uint32_t osbufferSize;
    
}*osStream;

typedef struct stream_stats_t{
    
    int32_t *alphabet;
    uint32_t *counts;
    int32_t *alphaMap;
    uint8_t *alphaExist;
    uint32_t alphabetCard;
    uint32_t step;
    uint32_t n;
    
}*stream_stats;

typedef struct arithStream_t{
    stream_stats* stats;
    Arithmetic_code a;
    osStream os;
}*arithStream;

typedef struct fasta_compressor_t{
    arithStream Ints;
    arithStream Chars;
    arithStream Signs;
}*fasta_compressor;





uint8_t read_bit_from_stream(osStream is);
uint32_t read_uint32_from_stream(uint32_t numBits, osStream is);
uint32_t send_bit_to_os(uint8_t bit, osStream os);
uint32_t send_uint32_to_os(uint32_t num, uint8_t numBits, osStream os);
osStream initialize_osStream(uint32_t osBuffer, FILE *fos, FILE *fosA, uint8_t inStream);

Arithmetic_code initialize_arithmetic_encoder(uint32_t m);
uint32_t arithmetic_encoder_step(Arithmetic_code a, stream_stats stats, int32_t x, osStream os);
int encoder_last_step(Arithmetic_code a, osStream os);
uint32_t arithmetic_decoder_step(Arithmetic_code a, stream_stats stats, osStream is);
uint32_t decoder_last_step(Arithmetic_code a, stream_stats stats);

stream_stats* initialize_stream_stats_Ints();
stream_stats* initialize_stream_stats_Chars(uint8_t **BP_trans);
stream_stats* initialize_stream_stats_Signs();

uint32_t rescale_stats(stream_stats s);
uint32_t rescale_stats_var(stream_stats s);
uint32_t update_stats_Ints(stream_stats stats, int32_t x, uint32_t m, uint8_t stepSize);
uint32_t update_stats_Chars(stream_stats stats, int32_t x, uint32_t m);
uint32_t update_stats_Signs(stream_stats stats, int32_t x, uint32_t m);

uint8_t extract_byte(unsigned int x, uint8_t B);

uint32_t compress_Signs(arithStream I, uint32_t x);
uint32_t decompress_Signs(arithStream I);
uint32_t compress_Int_Byte(arithStream I, uint32_t x, unsigned int B);
uint32_t decompress_Int_Byte(arithStream I, unsigned int B);
void compress_Ints(arithStream I, uint32_t x);
uint32_t decompress_Ints(arithStream I);
uint32_t compress_Chars(arithStream I, char target, char ref);
char decompress_Chars(arithStream I, char ref);

#endif
