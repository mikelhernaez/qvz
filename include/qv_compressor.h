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

#define OS_STREAM_BUF_LEN		4096

#define COMPRESSION 0
#define DECOMPRESSION 1

typedef struct Arithmetic_code_t{
    
    uint32_t scale3;
    uint32_t m;
    uint32_t l;
    uint32_t u;
    uint32_t t;
}*Arithmetic_code;

typedef struct os_stream_t {
	FILE *fp;
	uint8_t *buf;
	uint32_t bufPos;
	uint8_t bitPos;
	uint64_t written;
} *osStream;

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




// Stream interface
struct os_stream_t *alloc_os_stream(FILE *fp, uint8_t in);
void free_os_stream(struct os_stream_t *);
uint8_t stream_read_bit(struct os_stream_t *);
uint32_t stream_read_bits(struct os_stream_t *os, uint8_t len);
void stream_write_bit(struct os_stream_t *, uint8_t);
void stream_write_bits(struct os_stream_t *os, uint32_t dw, uint8_t len);
void stream_finish_byte(struct os_stream_t *);
void stream_write_buffer(struct os_stream_t *);

// Arithmetic ncoder interface
Arithmetic_code initialize_arithmetic_encoder(uint32_t m);
uint32_t arithmetic_encoder_step(Arithmetic_code a, stream_stats stats, int32_t x, osStream os);
int encoder_last_step(Arithmetic_code a, osStream os);
uint32_t arithmetic_decoder_step(Arithmetic_code a, stream_stats stats, osStream is);
uint32_t decoder_last_step(Arithmetic_code a, stream_stats stats);

stream_stats** initialize_stream_stats(struct cond_quantizer_list_t *q_list);

uint32_t update_stats(stream_stats stats, int32_t x, uint32_t m);

void compress_qv(arithStream as, uint32_t x, uint32_t column, uint32_t idx);
uint32_t decompress_qv(arithStream as, uint32_t column, uint32_t idx);

qv_compressor initialize_qv_compressor(char osPath[], uint8_t streamDirection, struct cond_quantizer_list_t *q_list);

uint32_t start_qv_compression(FILE *fp, char* osPath, struct cond_quantizer_list_t *qlist, double *dis);
uint32_t start_qv_decompression(FILE *fop, char* isPath, struct cond_quantizer_list_t *qlist);

#endif
