//
//  os_stream.c
//  samSecComp
//
//  Created by Mikel Hernaez on 6/2/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "qv_compressor.h"




uint8_t read_bit_from_stream(osStream is){
    
    uint8_t bit;
    
    bit = is->bufferByte << is->bitPos++;
    bit >>= 7;
    
    if (is->bitPos == 8)
        is->bufferByte = fgetc(is->fos), is->bitPos = 0;
    
    return bit;
    
}
uint32_t read_uint32_from_stream(uint32_t numBits, osStream is){
    
    uint32_t num = 0;
    int32_t bitPos;
    uint8_t tempBit = 0;
    
    for (bitPos = numBits - 1; bitPos >= 0; bitPos--){
        tempBit = read_bit_from_stream(is);
        num |= tempBit << bitPos;
    }
    
    return num;
    
    
}
uint32_t send_bit_to_os(uint8_t bit, osStream os){
    
    bit <<= --os->bitPos;
    os->bufferByte |= bit;
    
    if (os->bitPos == 0) {
        fputc(os->bufferByte, os->fos), os->bitPos = 8, os->bufferByte = 0;
        os->osLength++;
    }
    
    return os->osLength;
    
}

uint32_t send_uint32_to_os(uint32_t num, uint8_t numBits, osStream os){
    
    uint32_t tmp = 0;
    int bitPos;
    uint8_t tempBit = 0;
    
    for (bitPos = numBits - 1; bitPos >= 0; bitPos--){
        tmp = num;
        tmp <<= 31 - bitPos;
        tmp >>= 31;
        tempBit = (uint8_t)tmp;
        send_bit_to_os(tempBit, os);
    }
    
    return os->osLength;
    
}

osStream initialize_osStream(uint32_t osBuffer, FILE *fos, FILE *fosA, uint8_t inStream){
    
    osStream os;
    
    os = (osStream) calloc(1, sizeof(struct osStream_t));
    
    os->fos = fos;
    os->osBuffer = (uint8_t*) calloc(osBuffer, sizeof(uint8_t));
    os->osbufferSize = osBuffer;
    os->bitPos = 8;
    
    if (inStream) {
        os->bitPos = 0;
        os->bufferByte = fgetc(fos);
    }
    
    return os;
}