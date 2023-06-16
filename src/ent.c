/* @file  ent.c
**
** @@
******************************************************************************/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <slow5/slow5.h>
#include "error.h"
#include "sigtk.h"

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
    {"no-header",no_argument,0,'n'},               //3 suppress header
};

double entropy(int16_t *raw_signal, uint64_t len_raw_signal){
    int64_t *counts = (int64_t *)malloc(65536*sizeof(int64_t));
    MALLOC_CHK(counts);
    memset(counts,0,65536*sizeof(int64_t));
    int64_t total = 0;

    for(uint64_t i=0;i<len_raw_signal;i++){
        uint16_t a= raw_signal[i];
        counts[a]++;
        //double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
        //printf("%f ",pA);
    }

    double ent = 0;
    for(uint64_t i=0;i<65536;i++){
        if(counts[i]>0){
            double p = (double)counts[i]/(double)len_raw_signal;
            ent -= p*log2(p);
        }
        total+=counts[i];
    }
    assert(total==len_raw_signal);
    free(counts);

    return ent;
}

static inline uint32_t _zigzag_encode_32 (int32_t val) {
	return (val + val) ^ (val >> 31);
}

void zigzag_delta_encode(const int32_t * in, uint32_t * out, size_t N, int32_t prev) {
    for (size_t i = 0; i < N; i++) {
      out[i] = _zigzag_encode_32(in[i] - prev);
      prev = in[i];
    }
}

int entmain(int argc, char* argv[]) {

    const char* optstring = "hV";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    int8_t hdr = 1;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigtk %s\n",SIGTK_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='n'){
            hdr = 0;
        }
    }

    if (argc-optind!=1 || fp_help == stdout) {
        fprintf(fp_help,"Usage: sigtk ent a.blow5\n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   -n                         suppress header\n");
        fprintf(fp_help,"   --version                  print version\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    slow5_file_t *sp = slow5_open(argv[optind],"r");
    if(sp==NULL){
       fprintf(stderr,"Error in opening file\n");
       exit(EXIT_FAILURE);
    }
    slow5_rec_t *rec = NULL;
    int ret=0;

    if(hdr) printf("read_id\traw_ent\tdelta_ent\tbyte_ent\n");

    while((ret = slow5_get_next(&rec,sp)) >= 0){
        printf("%s\t",rec->read_id);

        //raw
        double ent = entropy(rec->raw_signal,rec->len_raw_signal);
        printf("%f",ent);

        //zigzag delta
        int32_t *in = (int32_t *)malloc(rec->len_raw_signal*sizeof(int32_t));
        MALLOC_CHK(in);
        for(uint64_t i=0;i<rec->len_raw_signal;i++){
            in[i] = rec->raw_signal[i];
        }
        uint32_t *delta = (uint32_t *)malloc(rec->len_raw_signal*sizeof(uint32_t));
        MALLOC_CHK(delta);
        zigzag_delta_encode(in,delta,rec->len_raw_signal,0);

        int16_t *out = (int16_t *)malloc(rec->len_raw_signal*sizeof(int16_t));
        MALLOC_CHK(out);
        for(uint64_t i=0;i<rec->len_raw_signal;i++){
            out[i] = delta[i];
        }

        ent = entropy(out,rec->len_raw_signal-1);
        printf("\t%f",ent);

        //byte split
        int16_t *out2 = (int16_t *)malloc(rec->len_raw_signal*sizeof(int16_t));
        MALLOC_CHK(out2);
        int16_t *out3 = (int16_t *)malloc(rec->len_raw_signal*sizeof(int16_t));
        MALLOC_CHK(out3);

        for(uint64_t i=0;i<rec->len_raw_signal-1;i++){
            uint16_t a= out[i];
            assert(a/256 == a>>8);
            assert(a%256 == (uint16_t)(a<<8)>>8);

            out2[i] = a/256;
            out3[i] = a%256;
        }
        ent = entropy(out2,rec->len_raw_signal-1)+entropy(out3,rec->len_raw_signal-1);
        printf("\t%f",ent);

        for(uint64_t i=0;i<rec->len_raw_signal-1;i++){
            int32_t a= out2[i];
            int32_t b = out3[i];
            int32_t c = a*256+b;
            assert(c==out[i]);
        }

        free(in);
        free(delta);
        free(out);
        free(out2);
        free(out3);

        printf("\n");
    }

    if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
        fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
        exit(EXIT_FAILURE);
    }

    slow5_rec_free(rec);

    slow5_close(sp);

    return 0;

}
