//
//  fasta_compressor.h
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef qv_compressor_h

#define qv_compressor_h

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

#include "codebook.h"

#define m_arith  22


#define COMPRESSION 0
#define DECOMPRESSION 1

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
    
    uint32_t *counts;
    uint32_t alphabetCard;
    uint32_t step;
    uint32_t n;
    
}*stream_stats;

typedef struct arithStream_t{
    stream_stats** stats;
    Arithmetic_code a;
    osStream os;
}*arithStream;

typedef struct qv_compressor_t{
    arithStream Quals;
}*qv_compressor;





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

stream_stats** initialize_stream_stats(struct cond_quantizer_list_t *q_list);

uint32_t update_stats(stream_stats stats, int32_t x, uint32_t m);

uint32_t compress_qv(arithStream as, uint32_t x, uint32_t column, uint32_t idx);
uint32_t decompress_qv(arithStream as, uint32_t column, uint32_t idx);

qv_compressor initialize_qv_compressor(char osPath[], uint8_t streamDirection, struct cond_quantizer_list_t *q_list);

uint32_t start_qv_compression(FILE *fp, char* osPath, struct cond_quantizer_list_t *qlist, double *dis);
uint32_t start_qv_decompression(FILE *fop, char* isPath, struct cond_quantizer_list_t *qlist);

#endif
