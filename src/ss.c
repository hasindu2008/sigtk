/* @file  ss.c
**
** @@
******************************************************************************/

#define _XOPEN_SOURCE 700
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include "error.h"
#include "sigtk.h"

#define MAX_LEN_KMER 2000

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
};

typedef struct{
    char *rid;
    int32_t qlen;
    int32_t query_start;
    int32_t query_end;
    int8_t strand;
    char *tid;
    int32_t tlen;
    int32_t target_start;
    int32_t target_end;
    uint8_t mapq;
    char *ss;
}paf_rec_t;

paf_rec_t *parse_paf_rec(char *buffer){

    char *pch=NULL;

    paf_rec_t *paf = (paf_rec_t *)malloc(sizeof(paf_rec_t));
    MALLOC_CHK(paf);

    //read name
    pch = strtok (buffer,"\t\r\n"); assert(pch!=NULL);
    paf->rid = strdup(pch);

    //readlen
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->qlen = atoi(pch);

    //query start
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->query_start = atoi(pch);

    //query end
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->query_end= atoi(pch);

    //relative strand
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    if(strcmp(pch,"+")==0){
        paf->strand=0;
    }
    else if(strcmp(pch,"-")==0){
        paf->strand=1;
    }
    else{
        assert(0);
    }

    //targetname
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->tid = strdup(pch);

    //target len
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->tlen = atoi(pch);

    //target start
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->target_start = atoi(pch);

    //target end
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->target_end= atoi(pch);

    //num residue
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

    //num block
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

    //mapq
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->mapq = atoi(pch);

    paf->ss = NULL;
    while((pch = strtok(NULL,"\t\r\n"))){
        if(strncmp("ss:Z:",pch,5)==0){
            //fprintf(stderr,"ss:Z:%s\n",pch);
            paf->ss = strdup(pch+5);
        }
    }
    //error check
    if(paf->ss==NULL){
        ERROR("ss:Z: tag not found in paf record for %s",paf->rid);
        exit(EXIT_FAILURE);
    }

    return paf;
}

void free_paf_rec(paf_rec_t *paf){
    free(paf->rid);
    free(paf->tid);
    free(paf->ss);
    free(paf);
}

static void ss2tsv(paf_rec_t *paf){

    char *ss=paf->ss;
    int start_raw=paf->query_start; int end_raw=paf->query_end; //int len_raw_signal=paf->qlen;
    int start_kmer=paf->target_start; int end_kmer=paf->target_end; int len_kmer=paf->tlen;

    // Raw signal start index for the corresponding k-mer and Raw signal end index for the corresponding k-mer
    size_t cap = MAX_LEN_KMER;
    int *st_raw_idx = (int *)malloc(sizeof(int)*cap);
    MALLOC_CHK(st_raw_idx);
    int *end_raw_idx = (int *)malloc(sizeof(int)*cap);
    MALLOC_CHK(end_raw_idx);

    //intialise to -1
    for(int i=0; i<cap; i++){ st_raw_idx[i]=end_raw_idx[i]=-1; }

    int st_k = start_kmer; int end_k = end_kmer; //if DNA, start k-kmer index is start_kmer column in paf and end k-kmer index is end_kmer column in paf
    int8_t rna = start_kmer > end_kmer ? 1 : 0; //if RNA start_kmer>end_kmer in paf
    if(rna){ st_k = end_kmer; end_k = start_kmer; } //if RNA, start k-kmer index is end_kmer column in paf and end k-kmer index is start_kmer column in paf

    int i_k = st_k; int i_raw = start_raw; //current k-mer index and current raw signal index

    //buffer for storing digits preceding each operation and its index
    char buff[11]; int i_buff=0;

    while(*ss){

        if(*ss==',' || *ss=='I' || *ss=='D'){
            if(i_buff <= 0){ fprintf(stderr,"Bad ss: Preceding digit missing\n"); exit(1); }//if nothing in buff

            buff[i_buff]=0; //null terminate buff
            int num = atoi(buff);
            if(num < 0){ fprintf(stderr,"Bad ss: Cannot have negative numbers\n"); exit(1); }
            i_buff=0; buff[0]=0; //reset buff

            if(*ss=='I'){ //if an insertion, current raw signal index is incremented by num
                i_raw += num;
            } else if(*ss=='D'){ //if an deletion, current k-mer index is incremented by num
                i_k += num;
            } else if (*ss==','){ //if a mapping, increment accordingly and set raw signal indices for the current k-mer
                end_raw_idx[i_k] = i_raw; i_raw += num;
                st_raw_idx[i_k] = i_raw; i_k++;
            }
            if(i_k >= cap){
                st_raw_idx=(int *)realloc(st_raw_idx, sizeof(int)*cap*2);
                MALLOC_CHK(st_raw_idx);
                end_raw_idx=(int *)realloc(end_raw_idx, sizeof(int)*cap*2);
                MALLOC_CHK(end_raw_idx);
                for(int i=cap; i<cap*2; i++){ st_raw_idx[i]=end_raw_idx[i]=-1; }
                cap *= 2;
            }
        } else {
            if(!isdigit(*ss)){ fprintf(stderr,"Bad ss: A non-digit found when expected a digit\n"); exit(1); }
            buff[i_buff++]=*ss;
        }
        ss++;
    }

    if(i_raw!=end_raw){ fprintf(stderr,"Bad ss: Signal end mismatch\n"); exit(1); } //current raw signal index should be equal to end_raw
    if(i_k!=end_k){ fprintf(stderr,"Bad ss: Kmer end mismatch\n"); exit(1); } //current k-mer index should be equal to end_k

    for(int i=st_k; i<end_k; i++){
        if(end_raw_idx[i]==-1){
            if(st_raw_idx[i] != -1) { fprintf(stderr,"Bad ss: This shoud not have happened\n"); exit(1); }//if st_raw_idx[i] is -1, then end_raw_idx[i] should also be -1
            printf("%s\t%d\t.\t.\n", paf->rid, rna ? len_kmer-i-1 : i);
        }else {
            printf("%s\t%d\t%d\t%d\n", paf->rid, rna ? len_kmer-i-1 : i, end_raw_idx[i], st_raw_idx[i]);
        }
    }

    free(st_raw_idx);
    free(end_raw_idx);

}

static void do_paf2tsv(char *paffile){

    //buffers for getline
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char)*(bufferSize));
    MALLOC_CHK(buffer);
    int readlinebytes=1;

    FILE *fp = fopen(paffile,"r");
    F_CHK(fp,paffile);

    printf("read_id\tkmer_idx\tstart_raw_idx\tend_raw_idx\n");

	while(1){
        readlinebytes=getline(&buffer, &bufferSize, fp);
        if(readlinebytes == -1){
            break;
        }

        paf_rec_t *paf = parse_paf_rec(buffer);

        ss2tsv(paf);

        free_paf_rec(paf);
    }

    free(buffer);
    fclose(fp);

    return;

}


int ssmain(int argc, char* argv[]) {

    const char* optstring = "hV";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    // int8_t hdr = 1;

    // opt_t opt = {0,0,0};

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigtk %s\n",SIGTK_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        }
    }


    if (argc-optind!=2 ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigtk ss paf2tsv in.paf\n");
        //fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    if(strcmp(argv[optind],"paf2tsv") == 0){
        do_paf2tsv(argv[optind+1]);
    }

    return 0;

}
