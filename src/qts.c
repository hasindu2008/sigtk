/* @file  qts.c
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
    {"bits",required_argument,0,'b'},              //4 number of LSB bits
    {"method",required_argument,0,'m'}             //5 quantization method
};

int round_to_power_of_2(int number, int number_of_bits) {
    //create a binary mask with the specified number of bits
    int bit_mask = (1 << number_of_bits) - 1;

    //extract out the value of the LSBs considered
    int lsb_bits = number & bit_mask;


    int round_threshold = (1 << (number_of_bits - 1));

    //check if the least significant bits are closer to 0 or 2^n
    if (lsb_bits < round_threshold) {
        return (number & ~bit_mask) + 0; //round down to the nearest power of 2
    } else {
        return (number & ~bit_mask) + (1 << number_of_bits); //round up to the nearest power of 2
    }
}



int qtsmain(int argc, char* argv[]) {

    const char* optstring = "hVv:o:b:m:";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    char *out_fn = NULL;

    int8_t b = 1; //number of LSB bits to truncate

    char *method = "round"; //quantization method

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigtk %s\n",SIGTK_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='o'){
            out_fn = optarg;
        } else if (c=='b'){
            b = atoi(optarg);
            if (b < 1 || b > 8) {
                fprintf(stderr, "Error: number of bits to truncate must be between 1 and 8\n");
                exit(EXIT_FAILURE);
            }
        } else if (c=='m'){
            method = optarg;
        }
    }

    if (argc-optind!=1 || fp_help == stdout) {
        fprintf(fp_help,"Usage: sigtk qts a.blow5 -o out.blow5\n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                            help\n");
        fprintf(fp_help,"   -o FILE                       output file\n");
        fprintf(fp_help,"   --version                     print version\n");
        fprintf(fp_help,"   -b INT                        number of lower significant bits to eliminate [%d]\n", b);
        fprintf(fp_help,"   -m [floor|round|fill-ones]    quantisation method [round]\n");
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

    while((ret = slow5_get_next(&rec,sp)) >= 0){
        if(strcmp(method, "floor") == 0){// truncate b LSBs
            for(uint64_t i=0;i<rec->len_raw_signal;i++){
                rec->raw_signal[i] = ( rec->raw_signal[i] >> b ) << b;
            }
        } else if(strcmp(method, "round") == 0){
            for(uint64_t i=0;i<rec->len_raw_signal;i++){
                rec->raw_signal[i] = round_to_power_of_2(rec->raw_signal[i], b);
            }
        } else if(strcmp(method, "fill-ones") == 0) { //fill the interested LSBs with 1s
            for(uint64_t i=0;i<rec->len_raw_signal;i++){
                rec->raw_signal[i] = (rec->raw_signal[i]) | ((1 << b)-1);
            }
        } else {
            fprintf(stderr,"Unknown method for -m. Available options are floor,round,fill-ones.\n");
            exit(EXIT_FAILURE);
        }

        //write to file
        if(slow5_write(rec, sp_w) < 0){
            fprintf(stderr,"Error writing record!\n");
            exit(EXIT_FAILURE);
        }

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
