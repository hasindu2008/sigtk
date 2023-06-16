/* @file  misc.c
**
** @@
******************************************************************************/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "sigtk.h"
#include "misc.h"
#include "error.h"

float *signal_in_picoamps(slow5_rec_t *rec){
    int16_t* rawptr = rec->raw_signal;
    float range = rec->range;
    float digitisation = rec->digitisation;
    float offset = rec->offset;
    int32_t nsample = rec->len_raw_signal;

    // convert to pA
    float *current_signal = (float*)malloc(sizeof(float) * nsample);
    MALLOC_CHK(current_signal);

    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        current_signal[j] = ((float)rawptr[j] + offset) * raw_unit;
    }

    return current_signal;
}

int8_t drna_detect(slow5_file_t *sp){

    const slow5_hdr_t* hdr = sp->header;
    int8_t rna = 0;
    char *exp =slow5_hdr_get("experiment_type", 0, hdr);
    if(exp==NULL){
        WARNING("%s","experiment_type not found in SLOW5 header. Assuming genomic_dna");
        return 0;
    }
    if (strcmp(exp,"genomic_dna")==0){
        rna = 0;
    }else if (strcmp(exp,"rna")==0){
        rna = 1;
        INFO("%s","RNA data detected.");
    } else {
        WARNING("Unknown experiment type: %s. Assuming genomic_dna", exp);
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, exp)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Defaulted to %s", curr, exp, i, exp);
        }
    }
    return rna;
}

void drna_mismatch(slow5_file_t *sp, int8_t rna){
    const char *expected = rna ? "rna" : "genomic_dna";
    const slow5_hdr_t* hdr = sp->header;
    for(uint32_t  i=0; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, expected)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Double check for --rna.", curr, expected, i);
        }
    }
}
