//
//  fasta_stream.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "qv_compressor.h"

uint8_t** generate_char_transitions(chromosome chr, arithStream signs){
    
    uint8_t** BP_transitions;
    
    BASEPAIR refBP, targetBP;
    
    uint32_t i = 0, j = 0;
    
    BP_transitions = (uint8_t **)calloc(35, sizeof(uint8_t*));
    
    for (i = 0; i < 35; i++) {
        BP_transitions[i] = (uint8_t *)calloc(35, sizeof(uint8_t));
    }
    
    // Store the rest of the numbers related to the edit distance between each of the pair of strings
    while (chr != NULL) {
        
        // Check instructions
        for (i = 0; i < chr->numInst - 1; i++) {
            
            refBP = char2BP(chr->instructions[i].refChar);
            targetBP = char2BP(chr->instructions[i].targetChar);
            
            BP_transitions[refBP][targetBP] = 1;
        }
        
        // Check insertions
        for (i = 0; i < chr->numInse; i++) {
            
            refBP = char2BP(chr->insertions[i].refChar);
            targetBP = char2BP(chr->insertions[i].targetChar);
            
            BP_transitions[refBP][targetBP] = 1;
        }
        
        // Check substitutions
        for (i = 0; i < chr->numSubs; i++) {
            
            refBP = char2BP(chr->substitutions[i].refChar);
            targetBP = char2BP(chr->substitutions[i].targetChar);
            
            BP_transitions[refBP][targetBP] = 1;
        }
        
        chr = chr->next;
    }
    
    for (i = 0; i < 35; i++) {
        for (j = 0; j < 35; j++) {
            compress_Signs(signs, BP_transitions[i][j]);
        }
    }
    
    return BP_transitions;
}

uint8_t** regenerate_char_transitions(arithStream signs){
    
    uint8_t** BP_transitions;
    
    uint32_t i = 0, j = 0;
    
    // Allocate memory
    BP_transitions = (uint8_t **)calloc(35, sizeof(uint8_t*));
    
    for (i = 0; i < 35; i++) {
        BP_transitions[i] = (uint8_t *)calloc(35, sizeof(uint8_t));
    }
    
    // regenerate the transitions
    for (i = 0; i < 35; i++) {
        for (j = 0; j < 35; j++) {
            BP_transitions[i][j] = decompress_Signs(signs);
        }
    }
    
    return BP_transitions;
}


arithStream initialize_arithStream_Ints(char* osPath, uint8_t decompressor_flag){
    
    arithStream as;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_ints.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats_Ints();
    as->a = initialize_arithmetic_encoder(m_INTS);
    as->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return as;
    
}

arithStream initialize_arithStream_Signs(char* osPath, uint8_t decompressor_flag){
    
    arithStream as;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_signs.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats_Ints();
    as->a = initialize_arithmetic_encoder(m_SIGNS);
    as->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return as;
    
}

arithStream initialize_arithStream_Chars(char* osPath, uint8_t decompressor_flag, uint8_t **BP_trans){
    
    arithStream as;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_char.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats_Chars(BP_trans);
    as->a = initialize_arithmetic_encoder(m_CHARS);
    as->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return as;
    
}

fasta_compressor initialize_fasta_compressor(char osPath[], uint8_t streamDirection, chromosome chr){
    
    fasta_compressor s;
    uint8_t** BP_transitions;
    
    s = calloc(1, sizeof(struct fasta_compressor_t));
    
    s->Signs = initialize_arithStream_Signs(osPath, streamDirection);
    
    s->Ints = initialize_arithStream_Ints(osPath, streamDirection);
    
    if (streamDirection == COMPRESSION) {
        BP_transitions = generate_char_transitions(chr, s->Signs);
    }
    else
        BP_transitions = regenerate_char_transitions(s->Signs);
    
    
    s->Chars = initialize_arithStream_Chars(osPath, streamDirection, BP_transitions);
    
    
    return s;
}
