/* @file  events_main.c
**
** @@
******************************************************************************/
#include "sigtk.h"
#include "misc.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int compact = 0;

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
    {"output",required_argument, 0, 'o'},          //3 output to a file [stdout]
    {"rna",no_argument,0,0},                       //4 if RNA
    {"no-header",no_argument,0,'n'},               //5 suppress header
    {"compact",no_argument,0,'c'},                 //6 compact output
    {0, 0, 0, 0}};

void event_func(slow5_rec_t *rec, int8_t rna);
void stat_func(slow5_rec_t *rec, int8_t rna);
void seg_func(slow5_rec_t *rec, int8_t rna);
void pa_func(slow5_rec_t *rec, int8_t rna);
void stat_hdr();
void seg_hdr();
void event_hdr();
void pa_hdr();
void jnn_hdr();
void jnn_func(slow5_rec_t *rec, int8_t rna);

int cmain(int argc, char* argv[], char *mode) {

    const char* optstring = "o:hVnc";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    int8_t rna = 0;
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
        } else if (c=='c'){
            compact = 1;
        // } else if (c == 0 && longindex == 4){ //if RNA
        //     rna = 1;
        }
    }


    if (argc-optind<1 ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigtk %s reads.blow5 read_id1 read_id2 .. \n", mode);
        fprintf(fp_help,"       sigtk %s reads.blow5\n", mode);
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   -n                         supress header\n");
        fprintf(fp_help,"   -c                         compact output\n");
        fprintf(fp_help,"   --version                  print version\n");
        //fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    slow5_file_t *slow5file = slow5_open(argv[optind], "r");
    if (!slow5file) {
        ERROR("cannot open %s. \n", argv[optind]);
        exit(EXIT_FAILURE);
    }
    rna = drna_detect(slow5file);

    slow5_rec_t *rec = NULL;
    int ret=0;

    void (*func)(slow5_rec_t*, int8_t) = NULL;

    if (strcmp(mode,"event") == 0){
        if(hdr) event_hdr();
        func = event_func;
    }else if (strcmp(mode,"stat") == 0){
        if(hdr) stat_hdr();
        func = stat_func;
    }else if (strcmp(mode,"seg") == 0){
        if(hdr) seg_hdr();
        func = seg_func;
    }else if (strcmp(mode,"jnn") == 0){
        if(hdr) jnn_hdr();
        func = jnn_func;
    }else if (strcmp(mode,"pa") == 0){
        if(hdr) pa_hdr();
        func = pa_func;
    } else {
        ERROR("unknown mode %s\n", mode);
        exit(EXIT_FAILURE);
    }

    if(argc-optind == 1){
        while((ret = slow5_get_next(&rec,slow5file)) >= 0){
           func(rec,rna);
        }

        if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
            fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
            exit(EXIT_FAILURE);
        }
    }
    else{

        int ret_idx = slow5_idx_load(slow5file);
        if (ret_idx < 0) {
            ERROR("Error loading index file for %s\n", argv[optind]);
            exit(EXIT_FAILURE);
        }

        for(int i=optind+1 ; i<argc ;i++){
            char *read_id = argv[i];
            fprintf(stderr,"Read ID %s\n",read_id);
            int ret = slow5_get(read_id, &rec, slow5file);
            if(ret < 0){
                ERROR("%s","Error when fetching the read\n");
                exit(EXIT_FAILURE);
            }

            func(rec,rna);

        }

        slow5_idx_unload(slow5file);

    }

    slow5_rec_free(rec);
    slow5_close(slow5file);

    return 0;
}
