//
//  stats.c
//  samSecComp
//
//  Created by Mikel Hernaez on 6/2/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//


#include "qv_compressor.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  INITIALIZATION                                      //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////


stream_stats* initialize_stream_stats_Signs(){
    
    stream_stats *s;
    
    s = (stream_stats*) calloc(1, sizeof(stream_stats));
        
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
    // Allocate memory
    s[0]->counts = (uint32_t*) calloc(3, sizeof(uint32_t));
        
    // An extra for the cumcounts
    s[0]->counts += 1;
        
    s[0]->alphabetCard = 2;
    
    s[0]->counts[0] = 1;
    s[0]->counts[1] = 1;
    s[0]->n = 2;
    
        
    // STEP
    s[0]->step = 1;
    
    return s;
    
}

stream_stats* initialize_stream_stats_Ints(){
    
    stream_stats *s;
    
    unsigned int j = 0, i = 0;
    
    s = (stream_stats*) calloc(4, sizeof(stream_stats));
    
    for (j = 0; j < 4; j++) {
        
        s[j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[j]->counts = (uint32_t*) calloc(257, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[j]->counts += 1;
        
        s[j]->alphabetCard = 256;
        
        
        s[j]->n = 0;
        for (i = 0; i < 256; i++) {
            s[j]->counts[i] = 1;
            s[j]->n ++;
        }
        
        // STEP
        s[j]->step = 8;
    }

    return s;
    
}


stream_stats* initialize_stream_stats_Chars(uint8_t **BP_trans){
    
    stream_stats* s;
    
    s = (stream_stats*) calloc(35, sizeof(stream_stats));
    
    uint32_t i = 0, j;
    
    // Stats for bp_A,bp_C,bp_G,bp_T,bp_N,bp_U,bp_R,bp_Y,bp_K,bp_M,bp_S,bp_W,bp_B,bp_D,bp_H, bp_V, bp_X, bp_a,bp_c,bp_g,bp_t,bp_n,bp_u,bp_r,bp_y,bp_k,bp_m,bp_s,bp_w,bp_b,bp_d,bp_h, bp_v, bp_x
    
    for (j = 0; j < 35; j++) {
        
        s[j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[j]->counts = (uint32_t*) calloc(36, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[j]->counts += 1;
        
        s[j]->n = 0;
        
        // Check the possible transitions
        // i is the value of the target
        // j is the value of the reference
        
        for (i = 0; i < 35; i++) {
            s[j]->counts[i] = (BP_trans[j][i] == 1)? 1:0;
            s[j]->n += s[j]->counts[i];
        }
        
        s[j]->alphabetCard = 35;
        
        // STEP
        s[j]->step = 8;
    }
 
    // Transition between complements, more common
    //A
    s[0]->counts[1] += 8;
    s[0]->counts[2] += 8;
    s[0]->n += 16;
    //C
    s[1]->counts[0] += 8;
    s[1]->counts[3] += 8;
    s[1]->n += 16;
    //G
    s[2]->counts[0] += 8;
    s[2]->counts[3] += 8;
    s[2]->n += 16;
    //T
    s[3]->counts[1] += 8;
    s[3]->counts[2] += 8;
    s[3]->n += 16;
    
    return s;
    
}


//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  UPDATE                                              //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

uint32_t update_stats_Ints(stream_stats stats, int32_t x, uint32_t m, uint8_t stepSize){
    
    int32_t i = 0;
    
    // Update the statistics
    if (stepSize)
        stats->counts[x]+= stats->step, stats->n+= stats->step;
    else
        stats->counts[x]+= 1, stats->n+= 1;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = -1; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}


uint32_t update_stats_Chars(stream_stats stats, int32_t x, uint32_t m){
    
    int32_t i = 0;
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = 0; i < (int32_t) stats->alphabetCard; i++){
            if (stats->counts[i]) {
                stats->counts[i] >>= 1, stats->counts[i]++;
                stats->n += stats->counts[i];
            }
        }
    }
    
    return 1;
}

uint32_t update_stats_Signs(stream_stats stats, int32_t x, uint32_t m){
    
    int32_t i = 0;
    // Update the statistics
    stats->counts[x]+= stats->step, stats->n+= stats->step;
    
    // Rescale if necessary
    if (stats->n >= (1 << (m - 3))){
        
        stats->n = 0;
        for (i = 0; i < (int32_t) stats->alphabetCard; i++){
            stats->counts[i] >>= 1, stats->counts[i]++;
            stats->n += stats->counts[i];
        }
    }
    
    return 1;
}
