//
//  arith.c
//  samSecComp
//
//  Created by Mikel Hernaez on 6/2/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <stdio.h>
#include "os_stream.c"
#include "stats.c"

Arithmetic_code initialize_arithmetic_encoder(uint32_t m){
    
    Arithmetic_code a_code;
    
    a_code = (Arithmetic_code) calloc(1, sizeof(struct Arithmetic_code_t));
    
    a_code->m = m;
    a_code->l = 0;
    a_code->u = 1;
    a_code->u <<= (a_code->m);
    a_code->u--;
    
    return a_code;
    
}

uint32_t arithmetic_encoder_step(Arithmetic_code a, stream_stats stats, int32_t x, osStream os){
    
    uint64_t range = 0;
    uint8_t msbU = 0, msbL = 0, E1_E2 = 0, E3 = 0, smsbL = 0, smsbU = 0;
    uint32_t osSize = 0, cumCountX = 0, cumCountX_1 = 0;
    int32_t i;
    
    range = a->u - a->l + 1;
    
    for (i = -1; i <= x; i++)
        cumCountX += stats->counts[i];
    
    cumCountX_1 = cumCountX - stats->counts[i-1];
    
    a->u = a->l + (uint32_t)( (range * cumCountX) / stats->n ) - 1;
    a->l = a->l + (uint32_t)( (range * cumCountX_1 ) / stats->n );
    
    // Check the rescaling conditions
    msbL = a->l >> (a->m - 1);
    msbU = a->u >> (a->m - 1);
    
    E1_E2 = (msbL == msbU);
    
    // If E1 or E2 doen't hold, check E3
    if (!E1_E2) {
        smsbL = a->l >> (a->m - 2);
        smsbU = (a->u >> (a->m - 2)) ^ 2;
        E3 = (smsbL == 1 && smsbU == 0);
    }
    else E3 = 0;
    
    // While any of E conditions hold
    while (E1_E2 || E3) {
        
        // If E1 or E2 holds rescale and recheck E1 and E2
        if (E1_E2) {
            send_bit_to_os(msbL, os);
            a->l <<= 1, a->l ^= msbL << a->m;
            a->u <<= 1, a->u ^= msbL << a->m, a->u++;
            
            while (a->scale3 > 0) {
                send_bit_to_os(!msbL, os);
                a->scale3--;
            }
            
            msbL = a->l >> (a->m - 1);
            msbU = a->u >> (a->m - 1);
            
            E1_E2 = (msbL == msbU);
            
            if (!E1_E2) {
                smsbL = a->l >> (a->m - 2);
                smsbU = (a->u >> (a->m - 2)) ^ 2;
                E3 = (smsbL == 1 && smsbU == 0);
            }
            else E3 = 0;
        }
        
        // if condition E3 holds, increment scale3 and rescale.
        if (E3){
            a->scale3++;
            a->u <<= 1, a->u ^= (1 << a->m) ^ (1 << (a->m - 1)), a->u++;
            a->l <<= 1, a->l ^= (1 << (a->m - 1));
            
            msbL = a->l >> (a->m - 1);
            msbU = a->u >> (a->m - 1);
            
            E1_E2 = (msbL == msbU);
            
            if (!E1_E2) {
                smsbL = a->l >> (a->m - 2);
                smsbU = (a->u >> (a->m - 2)) ^ 2;
                E3 = (smsbL == 1 && smsbU == 0);
            }
            else E3 = 0;
        }
    }
    
    return osSize;
}

int encoder_last_step(Arithmetic_code a, osStream os){
    
    uint8_t msbL = 0;
    uint32_t osSize;
    
    msbL = a->l >> (a->m - 1);
    
    // Write the msb of the tag (l)
    send_bit_to_os(msbL, os);
    
    // write as many !msbL as scale3 left
    while (a->scale3 > 0) {
        send_bit_to_os(!msbL, os);
        a->scale3--;
    }
    
    // write the rest of the tag (l)
    osSize = send_uint32_to_os(a->l, a->m - 1, os);
    
    fclose(os->fos);
    
    return osSize;
    
}

uint32_t arithmetic_decoder_step(Arithmetic_code a, stream_stats stats, osStream is){
    
    uint64_t range = 0, tagGap = 0;
    
    int32_t k = -1, x = -1, i;
    
    uint32_t subRange = 0, cumCountX = 0, cumCountX_1 = 0, cumCount = 0;
    
    uint8_t msbU = 0, msbL = 0, msbT = 0, E1_E2 = 0, E3 = 0, smsbL = 0, smsbU = 0;
    
    range = a->u - a->l + 1;
    tagGap = a->t - a->l + 1;
    
    subRange = (uint32_t)( (tagGap * stats->n - 1) / range );
    
    while (subRange >= cumCount)
        cumCount += stats->counts[k++];
    
    x = --k;
    
    for (i = -1; i <= x; i++)
        cumCountX += stats->counts[i];
    
    cumCountX_1 = cumCountX - stats->counts[i-1];
    
    a->u = a->l + (uint32_t)( (range * cumCountX) / stats->n ) - 1;
    a->l = a->l + (uint32_t)( (range * cumCountX_1) / stats->n);
    
    // Check the rescaling conditions.
    msbL = a->l >> (a->m - 1);
    msbU = a->u >> (a->m - 1);
    msbT = a->t >> (a->m - 1);
    
    E1_E2 = (msbL == msbU);
    
    // If E1 or E2 doen't hold, check E3
    if (!E1_E2) {
        smsbL = a->l >> (a->m - 2);
        smsbU = (a->u >> (a->m - 2)) ^ 2;
        E3 = (smsbL == 1 && smsbU == 0);
    }
    else E3 = 0;
    
    // While any of E conditions hold
    while (E1_E2 || E3) {
        
        // If E1 or E2 holds rescale and recheck E1, E2 and E3
        if (E1_E2) {
            a->l <<= 1, a->l ^= msbL << a->m;
            a->u <<= 1, a->u ^= msbL << a->m, a->u++;
            
            a->t <<= 1, a->t^=msbT << a->m, a->t += read_bit_from_stream(is);
            
            msbL = a->l >> (a->m - 1);
            msbU = a->u >> (a->m - 1);
            msbT = a->t >> (a->m - 1);
            
            E1_E2 = (msbL == msbU);
            
            if (!E1_E2) {
                smsbL = a->l >> (a->m - 2);
                smsbU = (a->u >> (a->m - 2)) ^ 2;
                E3 = (smsbL == 1 && smsbU == 0);
            }
            else E3 = 0;
        }
        
        // if condition E3 holds, increment scale3 and rescale.
        if (E3){
            
            a->l <<= 1, a->l ^= msbL << a->m;
            a->u <<= 1, a->u ^= msbU << a->m, a->u++;
            
            a->t <<= 1, a->t^=msbT << a->m, a->t += read_bit_from_stream(is);
            
            a->u ^= (1 << (a->m - 1));
            a->l ^= (1 << (a->m - 1));
            
            a->t ^= (1 << (a->m - 1));
            
            msbL = a->l >> (a->m - 1);
            msbU = a->u >> (a->m - 1);
            msbT = a->t >> (a->m - 1);
            
            E1_E2 = (msbL == msbU);
            
            if (!E1_E2) {
                smsbL = a->l >> (a->m - 2);
                smsbU = (a->u >> (a->m - 2)) ^ 2;
                E3 = (smsbL == 1 && smsbU == 0);
            }
            else E3 = 0;
        }
    }
    
    return x;
    
}

uint32_t decoder_last_step(Arithmetic_code a, stream_stats stats){
    
    uint64_t range, tagGap, subRange;
    uint32_t k = 0, cumCount = 0, x;
    
    range = a->u - a->l + 1;
    tagGap = a->t - a->l + 1;
    
    subRange = (tagGap * stats->n - 1) / range;
    
    while (subRange >= cumCount)
        cumCount += stats->counts[k++];
    
    x = --k;
    
    return x;
    
}

