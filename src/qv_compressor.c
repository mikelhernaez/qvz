//
//  fasta_compressor.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "qv_compressor.h"

void compress_qv(arithStream as, uint32_t x, uint32_t column, uint32_t idx){
    // Send x to the arithmetic encoder
    arithmetic_encoder_step(as->a, as->stats[column][idx], x, as->os);
    
    // Update the statistics
    update_stats(as->stats[column][idx], x, as->a->m);
}

uint32_t decompress_qv(arithStream as, uint32_t column, uint32_t idx){
    uint32_t x;
    
    // Send x to the arithmetic encoder
    x = arithmetic_decoder_step(as->a, as->stats[column][idx], as->os);
    // Update the statistics
    update_stats(as->stats[column][idx], x, as->a->m);
    
    return x;
}

uint32_t start_qv_compression(FILE *fp, char* osPath, struct cond_quantizer_list_t *qlist, double *dis) {
    unsigned int osSize = 0;
    
    qv_compressor qvc;
    
	uint32_t s = 0, idx = 0, lineCtr = 0, q_state = 0;
	double distortion = 0.0;
	uint32_t error = 0;
    uint8_t qv = 0, prev_qv = 0;
    uint32_t columns = qlist->columns;
    struct quantizer_t *q;

	char *line = (char *) _alloca(columns+3);
    
#ifdef DEBUG
    FILE *fref = fopen("fref.txt", "wt");
#endif
    
    // Initialize the compressor
    qvc = initialize_qv_compressor(osPath, COMPRESSION, qlist);
    
    // Start compressing the file
	distortion = 0.0;
	fgets(line, columns+2, fp);

	do {
        if(qlist->options->verbose && lineCtr%1000000 == 0){
            printf("Line: %dM\n", lineCtr/1000000);
        }
        lineCtr++;
        
		// Select first column's codebook with no left context
		q = choose_quantizer(qlist, 0, 0, &idx);
        
		// Quantize, compress and calculate error simultaneously
		// Note that in this version the quantizer outputs are 0-41, so the +33 offset is different from before
		qv = q->q[line[0]-33];
        q_state = get_symbol_index(q->output_alphabet, qv);
        compress_qv(qvc->Quals, q_state, 0, idx);
        
#ifdef DEBUG
		fputc(qv+33, fref);
#endif
        
		error = (line[0] - qv - 33)*(line[0] -qv - 33);
        prev_qv = qv;
        
		for (s = 1; s < columns; ++s) {
			// Quantize and compute error for MSE
			q = choose_quantizer(qlist, s, prev_qv, &idx);
			qv = q->q[line[s]-33];
            q_state = get_symbol_index(q->output_alphabet, qv);
            
#ifdef DEBUG
			fputc(qv+33, fref);
#endif
            
            compress_qv(qvc->Quals, q_state, s, idx);
			error += (line[s] - qv - 33)*(line[s] - qv - 33);
            prev_qv = qv;
		}
        
#ifdef DEBUG
       	fputc('\n', fref);
#endif
        
        distortion += error / ((double) columns);
        
		// Get next line from file
		fgets(line, columns+2, fp);
	} while (!feof(fp));
    
    osSize = encoder_last_step(qvc->Quals->a, qvc->Quals->os);
    
	qlist->lines = lineCtr;
    *dis = distortion / ((double) lineCtr);
    
    return osSize;
}

uint32_t start_qv_decompression(FILE *fop, char* isPath, struct cond_quantizer_list_t *qlist) {
    qv_compressor qvc;
    
	uint32_t s = 0, idx = 0, lineCtr = 0, q_state = 0;
    uint8_t prev_qv = 0;
    
    uint32_t columns = qlist->columns;
	uint32_t lines = qlist->lines;
    struct quantizer_t *q;

	char *line = (char *) _alloca(columns+2);
    line[columns] = '\n';
	line[columns+1] = '\0';
    
    // Initialize the compressor
    qvc = initialize_qv_compressor(isPath, DECOMPRESSION, qlist);
    
	// Last line has to be handled separately to clear the arithmetic decoder
	while (lineCtr < lines - 1) {
        if (qlist->options->verbose && lineCtr%1000000 == 0){
            printf("Line: %dM\n", lineCtr/1000000);
        }
        lineCtr++;
        
		// Select first column's codebook with no left context
		q = choose_quantizer(qlist, 0, 0, &idx);
        
		// Quantize, compress and calculate error simultaneously
		// Note that in this version the quantizer outputs are 0-41, so the +33 offset is different from before

        q_state = decompress_qv(qvc->Quals, 0, idx);
        line[0] = q->output_alphabet->symbols[q_state] + 33;
        prev_qv = line[0] - 33;
        
		for (s = 1; s < columns; ++s) {
			// Quantize and compute error for MSE
			q = choose_quantizer(qlist, s, prev_qv, &idx);
            q_state = decompress_qv(qvc->Quals, s, idx);
            line[s] = q->output_alphabet->symbols[q_state] + 33;
            prev_qv = line[s] - 33;
            
		}
        
        // Write this line to the output file, note '\n' at the end of the line buffer to get the right length
		fwrite(line, columns+1, sizeof(uint8_t), fop);
        
        // Write this line to the output file, note '\n' at the end of the line buffer to get the right length
		//for(int lala = 0; lala < 36; lala++)
        //    fprintf(fop, "%d ", line[lala]);
        //fputc('\n', fop);
        
	}
    
    // Last Line
    if(qlist->options->verbose && lineCtr%1000000 == 0){
        printf("Line: %dM\n", lineCtr/1000000);
    }
    lineCtr++;
    
    // Select first column's codebook with no left context
    q = choose_quantizer(qlist, 0, 0, &idx);
    
    // Quantize, compress and calculate error simultaneously
    // Note that in this version the quantizer outputs are 0-41, so the +33 offset is different from before
    
    q_state = decompress_qv(qvc->Quals, 0, 0);
    line[0] = q->output_alphabet->symbols[q_state] + 33;
    prev_qv = line[0] - 33;
    
    for (s = 1; s < columns - 1; ++s) {
        // Quantize and compute error for MSE
        q = choose_quantizer(qlist, s, prev_qv, &idx);
        q_state = decompress_qv(qvc->Quals, s, idx);
        line[s] = q->output_alphabet->symbols[q_state] + 33;
        prev_qv = line[s] - 33;
    }
    
    // Last column
    q = choose_quantizer(qlist, s, prev_qv, &idx);
    q_state = decoder_last_step(qvc->Quals->a, qvc->Quals->stats[s][idx]);
    line[s] = q->output_alphabet->symbols[q_state] + 33;
    
    // Write this line to the output file, note '\n' at the end of the line buffer to get the right length
    fwrite(line, columns+1, sizeof(uint8_t), fop);

	qlist->lines = lineCtr;
    
    return 0;
}
