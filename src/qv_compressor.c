//
//  fasta_compressor.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "qv_compressor.h"


uint32_t compress_Signs(arithStream I, uint32_t x){
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    // Send x to the arithmetic encoder
    arithmetic_encoder_step(I->a, I->stats[0], x, I->os);
    // Update the statistics
    update_stats_Signs(I->stats[0], x, I->a->m);
    
    return x;
}

uint32_t decompress_Signs(arithStream I){
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    uint32_t x = 0;
    
    // Send x to the arithmetic encoder
    x = arithmetic_decoder_step(I->a, I->stats[0], I->os);
    // Update the statistics
    update_stats_Signs(I->stats[0], x, I->a->m);
    
    return x;
}

uint32_t compress_Int_Byte(arithStream I, uint32_t x, unsigned int B){
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    // Send x to the arithmetic encoder
    arithmetic_encoder_step(I->a, I->stats[B], x, I->os);
    // Update the statistics
    update_stats_Ints(I->stats[B], x, I->a->m, BIG_STEP);
    
    return x;
}

uint32_t decompress_Int_Byte(arithStream I, unsigned int B){
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    unsigned int x = 0;
    
    // Send x to the arithmetic encoder
    x = arithmetic_decoder_step(I->a, I->stats[B], I->os);
    // Update the statistics
    update_stats_Ints(I->stats[B], x, I->a->m, BIG_STEP);
    
    return x;
}

uint8_t extract_byte(unsigned int x, uint8_t B){
    switch (B) {
        case 0:
            return (x & 0x000000ff);
        case 1:
            return (x & 0x0000ff00) >> 8;
        case 2:
            return (x & 0x00ff0000) >> 16;
        case 3:
            return (x & 0xff000000) >> 24;
            
        default:
            return 0;
    }
}

void compress_Ints(arithStream I, uint32_t x){
    
    compress_Int_Byte(I, extract_byte(x, 0), 0);
    compress_Int_Byte(I, extract_byte(x, 1), 1);
    compress_Int_Byte(I, extract_byte(x, 2), 2);
    compress_Int_Byte(I, extract_byte(x, 3), 3);
}

uint32_t decompress_Ints(arithStream I){
    
    uint8_t B0 = 0, B1 = 0, B2 = 0, B3 = 0;
    
    B0 = decompress_Int_Byte(I, 0);
    B1 = decompress_Int_Byte(I, 1);
    B2 = decompress_Int_Byte(I, 2);
    B3 = decompress_Int_Byte(I, 3);
    
    return ( (B3 << 24) ^ (B2 << 16) ^ (B1 << 8) ^ B0);
}

uint32_t compress_Chars(arithStream I, char target, char ref){
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    BASEPAIR Tbp, Rbp;
    
    Tbp = char2BP(target);
    Rbp = char2BP(ref);
    
    
    // Send x to the arithmetic encoder
    arithmetic_encoder_step(I->a, I->stats[Rbp], Tbp, I->os);
    // Update the statistics
    update_stats_Chars(I->stats[Rbp], Tbp, I->a->m);
    
    return 1;
}

char decompress_Chars(arithStream I, char ref){
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    BASEPAIR Tbp, Rbp;
    
    Rbp = char2BP(ref);
    
    // Send x to the arithmetic encoder
    Tbp = arithmetic_decoder_step(I->a, I->stats[Rbp], I->os);
    // Update the statistics
    update_stats_Chars(I->stats[Rbp], Tbp, I->a->m);
    
    return BP2char(Tbp);
}

uint32_t start_fasta_compression(chromosome chr, char* osPath, unsigned int numChr, unsigned int bpPerLine){
    
    unsigned int i = 0, prevInt = 0, osSize = 0;
    
    fasta_compressor fc;
    chromosome chrRoot;
    
    
    // Initialize the compressor
    fc = initialize_fasta_compressor(osPath, COMPRESSION, chr);
    
    // Start compressing the Ints
    
    // number of Chromosomes
    compress_Ints(fc->Ints, numChr);
    
    // number of base pairs per line
    compress_Ints(fc->Ints, bpPerLine);
    
    // [numInst numSubs numInse] x numChr
    chrRoot = chr;
    while (chr != NULL) {
        
        compress_Ints(fc->Ints, chr->numInst);
        compress_Ints(fc->Ints, chr->numInse);
        compress_Ints(fc->Ints, chr->numSubs);
        
        chr = chr->next;
    }
    chr = chrRoot;
    
    // Store the rest of the numbers related to the edit distance between each of the pair of strings
    while (chr != NULL) {
        
        // Compress instructions
        prevInt = 0;
        for (i = 0; i < chr->numInst - 1; i++) {
            
            compress_Ints(fc->Ints, abs(chr->instructions[i].pos - prevInt));
            
            compress_Ints(fc->Ints, chr->instructions[i].length);
            
            compress_Chars(fc->Chars, chr->instructions[i].targetChar, chr->instructions[i].refChar);

            // compress the signs
            compress_Signs(fc->Signs, (chr->instructions[i].pos < prevInt) );
            
            prevInt = chr->instructions[i].pos + chr->instructions[i].length;
        }
        
        // Compress the last Instruction without the chars (is always N to N)
        compress_Ints(fc->Ints, abs(chr->instructions[i].pos - prevInt));
        compress_Ints(fc->Ints, chr->instructions[i].length);
        compress_Signs(fc->Signs, (chr->instructions[i].pos < prevInt) );
        
        // Compress insertions
        prevInt = 0;
        for (i = 0; i < chr->numInse; i++) {
            
            compress_Ints(fc->Ints, abs(chr->insertions[i].pos - prevInt));
            
            compress_Chars(fc->Chars, chr->insertions[i].targetChar, chr->insertions[i].refChar);
            
            prevInt = chr->insertions[i].pos;
        }
        
        // Compress substitutions
        prevInt = 0;
        for (i = 0; i < chr->numSubs; i++) {
            
            compress_Ints(fc->Ints, abs(chr->substitutions[i].pos - prevInt));
            
            compress_Chars(fc->Chars, chr->substitutions[i].targetChar, chr->substitutions[i].refChar);
            
            prevInt = chr->substitutions[i].pos;
        }
        
        chr = chr->next;
    }

    osSize += encoder_last_step(fc->Ints->a, fc->Ints->os);
    osSize += encoder_last_step(fc->Signs->a, fc->Signs->os);
    osSize += encoder_last_step(fc->Chars->a, fc->Chars->os);
    
    return osSize;
    
}



