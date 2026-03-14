/* @file  filt.c
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
    {"output",required_argument,0,'o'},            //3 output file
    {"invert", no_argument, 0, 'i'},               //4 invert filter
    {"dist-thresh", required_argument, 0, 'd'},    //5 distance threshold
    {"trans-thresh", required_argument, 0, 't'},   //6 transition level


};

static int isto_write_read(slow5_rec_t *rec, slow5_file_t *sp, uint8_t invert, int dist_thresh, int trans_thresh){

    int high_count=0;

    int16_t val_prev=rec->raw_signal[0];
    for(int64_t i=1; i<rec->len_raw_signal; i++){
        int16_t val=rec->raw_signal[i];
        int16_t delta = val - val_prev;
        val_prev=val;
        int16_t abs = delta>=0 ? delta : -delta;
        if(abs>trans_thresh){
            high_count++;
        }
    }

    float gap = (float)rec->len_raw_signal/high_count;
    if (gap < dist_thresh ) {
        return (!invert);
    } else {
        return invert;
    }

}

int filtmain(int argc, char* argv[]) {

    const char* optstring = "hVv:o:id:t:";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    char *out_fn = NULL;

    uint8_t invert = 0;
    int dist_thresh = 100;
    int trans_thresh = 100;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigtk %s\n",SIGTK_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='o'){
            out_fn = optarg;
        } else if (c=='i'){
            invert = 1;
        } else if (c=='d'){
            dist_thresh = atoi(optarg);
        } else if (c=='t'){
            trans_thresh = atoi(optarg);
        }
    }

    if (argc-optind!=1 || fp_help == stdout) {
        fprintf(fp_help,"Usage: sigtk filt a.blow5 -o out.blow5\n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                            help\n");
        fprintf(fp_help,"   -o FILE                       output file\n");
        fprintf(fp_help,"   --invert                      invert filter\n");
        fprintf(fp_help,"   --dist-thresh VALUE           distance threshold (increase to be conservative when labelling reads as junk) [%d]\n", dist_thresh);
        fprintf(fp_help,"   --trans-thresh VALUE          transition level threshold [%d]\n", trans_thresh);
        fprintf(fp_help,"   --version                     print version\n");

        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    if(out_fn == NULL){
        fprintf(stderr,"Error: output file not specified\n");
        exit(EXIT_FAILURE);
    }

    //open the SLOW5 file for reading
    slow5_file_t *sp = slow5_open(argv[optind],"r");
    if(sp==NULL){
       fprintf(stderr,"Error in opening file\n");
       exit(EXIT_FAILURE);
    }

    //open the SLOW5 file for writing
    slow5_file_t *sp_w = slow5_open(out_fn, "w");
    if(sp_w==NULL){
        fprintf(stderr,"Error opening file!\n");
        exit(EXIT_FAILURE);
    }

    //backup the pointer to header in sp_w
    slow5_hdr_t *header=sp_w->header;
    sp_w->header = sp->header; //point the header from the input file to the output file

    if(slow5_hdr_write(sp_w) < 0){
        fprintf(stderr,"Error writing header!\n");
        exit(EXIT_FAILURE);
    }

    slow5_rec_t *rec = NULL;
    int ret=0;
    int written=0;
    int total = 0;

    while((ret = slow5_get_next(&rec,sp)) >= 0){
        total++;

        int eval = isto_write_read(rec, sp, invert, dist_thresh, trans_thresh);

        if(eval){
            //write to file
            if(slow5_write(rec, sp_w) < 0){
                fprintf(stderr,"Error writing record!\n");
                exit(EXIT_FAILURE);
            }
            written++;
        }

    }

    if(written==0){
        WARNING("%s\n","No reads were written. All were filtered out.");
    } else {
        fprintf(stderr, "%d/%d reads were written.\n", written, total);
    }

    if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
        fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
        exit(EXIT_FAILURE);
    }

    slow5_rec_free(rec);

    slow5_close(sp);

    //restore the header pointer
    sp_w->header = header;
    slow5_close(sp_w);


    return 0;

}