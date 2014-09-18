//
//  qv_stream.c
//  qvz
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "qv_compressor.h"


//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                            STATS UPDATE                                              //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

uint32_t update_stats(stream_stats stats, int32_t x, uint32_t m){
    
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


//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  INITIALIZATION                                      //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

stream_stats** initialize_stream_stats(struct cond_quantizer_list_t *q_list){
    
    stream_stats** s;
    
    s = (stream_stats**) calloc(q_list->columns, sizeof(stream_stats *));
    
    uint32_t i = 0, j = 0, k = 0;
    
    // Stats for all the different quantizers
    for (i = 0; i < q_list->columns; i++) {
        
        // Initialize stats for all the quantizers in column i (low and hi)
        s[i] = (stream_stats*) calloc(2*(q_list->input_alphabets[i]->size), sizeof(stream_stats));
        
        for (j = 0; j < 2*(q_list->input_alphabets[i]->size); j++){
            
            s[i][j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
            
            // Allocate memory for the counts
            s[i][j]->counts = (uint32_t*) calloc( (q_list->q[i][j]->output_alphabet->size) + 2, sizeof(uint32_t));
            
            // An extra for the cumcounts
            s[i][j]->counts += 1;
            
            s[i][j]->n = 0;
            
            // Initialize the quantizer's stats uniformly
            
            for (k = 0; k < q_list->q[i][j]->output_alphabet->size; k++) {
                s[i][j]->counts[k] = 1;
                s[i][j]->n++;
            }
            
            s[i][j]->alphabetCard = q_list->q[i][j]->output_alphabet->size;
            
            // STEP
            s[i][j]->step = 8;
            
        }
        
        
    }
    
    return s;
    
}


arithStream initialize_arithStream(char* osPath, uint8_t decompressor_flag, struct cond_quantizer_list_t *q_list){
    
    arithStream as;
    FILE *fp;
	uint32_t i;
    
    fp = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    if (decompressor_flag) {
        
        // First, read in the WELL state and set up the PRNG
        fread(q_list->well.state, sizeof(uint32_t), 32, fp);
    }
    
    else {
    
        // Initialize WELL state vector with libc rand (this initial vector needs to be copied to the decoder)
        srand((uint32_t) time(0));
        for (i = 0; i < 32; ++i) {
#ifndef DEBUG
            q_list->well.state[i] = rand();
#else
            qlist->well.state[s] = 0x55555555;
#endif
        }
        
        // Write the initial WELL state vector to the file first (fixed size of 32 bytes)
		// @todo strictly this needs to be stored in network order because we're interpreting it as a 32 bit int
		// but I am a bit too lazy for that right now
        fwrite(q_list->well.state, sizeof(uint32_t), 32, fp);
	}
	
	
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats(q_list);
    as->a = initialize_arithmetic_encoder(m_arith);
    as->os = initialize_osStream(1, fp, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    return as;
    
}

qv_compressor initialize_qv_compressor(char osPath[], uint8_t streamDirection, struct cond_quantizer_list_t *q_list){
    
    qv_compressor s;
    
    s = calloc(1, sizeof(struct qv_compressor_t));
    
    s->Quals = initialize_arithStream(osPath, streamDirection, q_list);
    
    
    return s;
}
