
/* @file genref.c
**

** @@
******************************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <getopt.h>

#include "sigtk.h"
#include "error.h"
#include "ref.h"
#include "misc.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

//#define NORMALISE_ALL 1

uint32_t read_model(model_t* model, const char* file, uint32_t type);
uint32_t set_model(model_t* model, uint32_t model_id);

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
    {"output",required_argument, 0, 'o'},          //3 output to a file [stdout]
    {"rna",no_argument,0,0},                       //4 if RNA
    {"no-header",no_argument,0,'n'},               //5 suppress header
    {0, 0, 0, 0}};

static inline void normalise(float *rawptr, uint64_t n){

    uint64_t start_idx =  0;
    uint64_t end_idx = n;

    float event_mean = 0;
    float event_var = 0;
    float event_stdv = 0;
    float num_samples = end_idx-start_idx;

    for(uint64_t j=start_idx; j<end_idx; j++){
        event_mean += rawptr[j];
    }
    event_mean /= num_samples;
    for(uint64_t j=start_idx; j<end_idx; j++){
        event_var += (rawptr[j]-event_mean)*(rawptr[j]-event_mean);
    }
    event_var /= num_samples;
    event_stdv = sqrt(event_var);

    for(uint64_t j=start_idx; j<end_idx; j++){
        rawptr[j] = (rawptr[j]-event_mean)/event_stdv;
    }

}

static inline void normalise_ref(refsynth_t *ref, int8_t rna){

    uint64_t n = 0;
    float mean = 0;
    float var = 0;
    float stdv = 0;
    for(int i=0; i<ref->num_ref; i++){
        n += ref->ref_lengths[i];
        for(int32_t j= 0; j<ref->ref_lengths[i]; j++){
            mean += ref->forward[i][j];
        }
    }
    mean /= n;

    for(int i=0; i<ref->num_ref; i++){
        for(int32_t j= 0; j<ref->ref_lengths[i]; j++){
            float a = ref->forward[i][j]-mean;
            var += a*a;
        }
    }
    var /= n;
    stdv =  sqrt(var);

    for(int i=0; i<ref->num_ref; i++){
        for(int32_t j= 0; j<ref->ref_lengths[i]; j++){
            ref->forward[i][j] = (ref->forward[i][j]-mean)/stdv;
            if (!rna){
                ref->reverse[i][j] = (ref->reverse[i][j]-mean)/stdv;
            }
        }

    }


}


void  gen_ref(const char *genome, model_t *pore_model, uint32_t kmer_size, uint32_t flag, int32_t query_size){

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(genome, "r");
    F_CHK(fp,genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);
    int8_t rna = flag & SIGTK_RNA;

    float *forward = NULL;
    float *reverse = NULL;

    int i = 0;
    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));

        int32_t ref_len;
        if(!rna || flag & SIGTK_REF){ //dna or use full reference
            ref_len = l+1-kmer_size;
        }
        else{ //rna
            uint32_t rlen_heu=query_size*1.5;
            ref_len = (rlen_heu > l+1-kmer_size ? l+1-kmer_size : rlen_heu);
            LOG_TRACE("Only %d bases of %d bases in reference sequence will be used\n", ref_len, l);
        }
        //int32_t ref_len = rna ? ((uint32_t)query_size > l+1-kmer_size ? l+1-kmer_size : query_size) : l+1-kmer_size;

        //int32_t ref_len =  l+1-kmer_size;

        forward = (float *) malloc(ref_len*sizeof(float));

        char *rc = NULL;
        if(!rna){
            reverse = (float *) malloc(ref_len*sizeof(float));
            rc = reverse_complement(seq->seq.s);
        }

        // fprintf(stderr,"%s\n",seq->seq.s);
        // fprintf(stderr,"%s\n",rc);

        if(!rna){ //dna (if DNA we just use the full reference)
            for (int j=0; j< ref_len; j++){
                uint32_t kmer_rank = get_kmer_rank(seq->seq.s+j, kmer_size);
                forward[j] = pore_model[kmer_rank].level_mean;

                kmer_rank = get_kmer_rank(rc+j, kmer_size);
                reverse[j] = pore_model[kmer_rank].level_mean;
            }
        }else{ //rna
            if(flag & SIGTK_INV){ //would be incorrect - have not tested recently
                VERBOSE("%s","Reversing the reference to be 5' -> 3'\n");
                // char *f = seq->seq.s;
                // char *reverse = (char *) malloc(strlen(f)*sizeof(char));
                // for(int j=0; j<(int)strlen(f); j++){
                //     reverse[j] = f[strlen(f)-j-1];
                // }
                char *seq_end = seq->seq.s+l-ref_len-(kmer_size-1);
                for(int j=0; j<ref_len; j++){
                    uint32_t kmer_rank = get_kmer_rank(seq_end+j, kmer_size);
                    forward[ref_len-j-1] = pore_model[kmer_rank].level_mean;
                }
                // for(int j=0; j<ref_len; j++){
                //     uint32_t kmer_rank = get_kmer_rank(seq->seq.s+j, kmer_size);
                //     ref->forward[i][ref_len-j-1] = pore_model[kmer_rank].level_mean;
                // }
            }
                // free(reverse);
            else{
                char *st;
                if(flag &  SIGTK_END){ //if the end of query then it is the beginning of the reference in RNA
                    st = seq->seq.s;

                }
                else{ //if the beginning of query then it is the end of the reference in RNA
                    st = seq->seq.s+l-ref_len-(kmer_size-1);
                }
                LOG_TRACE("%s:%ld-%ld\n",seq->name.s,st-seq->seq.s,st-seq->seq.s+ref_len);
                for (int j=0; j< ref_len; j++){
                    uint32_t kmer_rank = get_kmer_rank(st+j, kmer_size);
                    forward[j] = pore_model[kmer_rank].level_mean;
                }
            }
        }

        printf("%s\t%d\t+\t%d\t",seq->name.s,l,ref_len);
        for (int j=0; j< ref_len; j++){
            if (j < ref_len-1) printf("%f,",forward[j]); else printf("%f\n",forward[j]);
        }
        free(forward);
        if(!rna){
            printf("%s\t%d\t-\t%d\t",seq->name.s,l,ref_len);
            for (int j=0; j< ref_len; j++){
                if (j < ref_len-1) printf("%f,",reverse[j]); else printf("%f\n",reverse[j]);
            }

            free(reverse);
            free(rc);
        }

        i++;

    }

    kseq_destroy(seq);
    gzclose(fp);



}


int srefmain(int argc, char* argv[]){

    const char* optstring = "o:hVn";

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
        } else if (c == 0 && longindex == 4){ //if RNA
            rna = 1;
        }
    }


    if (argc-optind<1 ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigtk sref ref.fa \n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   -n                         suppress header\n");
        fprintf(fp_help,"   --version                  print version\n");
        fprintf(fp_help,"   --rna                      use RNA model\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    uint32_t flag = 0;
    if(rna){
        flag |= SIGTK_RNA;
        flag |= SIGTK_REF;
    }

    model_t *model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER);
    MALLOC_CHK(model);

    uint32_t kmer_size=0;
    if(flag & SIGTK_RNA){
        INFO("%s","builtin RNA nucleotide model loaded");
        kmer_size=set_model(model, MODEL_ID_RNA_NUCLEOTIDE);
    }
    else{
        kmer_size=set_model(model, MODEL_ID_DNA_NUCLEOTIDE);
    }


    if(hdr){
        printf("ref_name\tref_len\tstrand\tsig_len\tsig_mean\n");
    }

    gen_ref(argv[optind], model,kmer_size,flag, 250);

    free(model);

    return 0;

}