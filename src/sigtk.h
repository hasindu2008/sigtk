/* @file sigtk.h
**
******************************************************************************/

#ifndef SIGTK_H
#define SIGTK_H

#include <stdint.h>
#include "slow5/slow5.h"

#define SIGTK_VERSION "0.1.0"

//model types
#define MODEL_TYPE_NUCLEOTIDE 1
#define MODEL_TYPE_METH 2

#define MAX_KMER_SIZE 6 //maximum k-mer size
#define MAX_NUM_KMER 4096   //maximum number of k-mers in nucleotide model
#define MAX_NUM_KMER_METH 15625 //maximum number of k-mers in methylated model

//default model IDs
#define MODEL_ID_DNA_NUCLEOTIDE 1
#define MODEL_ID_RNA_NUCLEOTIDE 2

/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/

#define SIGTK_RNA 0x001 //if RNA or not
#define SIGTK_DTW 0x002 //if unused
#define SIGTK_INV 0x004 //if set, reverse reference events instead of query events
#define SIGTK_SEC 0x008 //if secondaries are printed
#define SIGTK_REF 0x010 //map to the whole reference
#define SIGTK_END 0x020 //map the end of the query

#define SECONDARY_CAP 5 //maximum number of secondary events to print

#define WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define STEAL_THRESH 1 //stealing threshold

//linear segment alignment record
#define LSAR_TJUMP 'J'      //jump in the target
#define LSAR_QJUMP 'j'      //jump in the query
#define LSAR_TMATCH 'M'      //match in the target
#define LSAR_QMATCH 'm'   //match in the query

/* a single signal-space event : adapted from taken from scrappie */
typedef struct {
    uint64_t start;
    float length; //todo : cant be made int?
    float mean;
    float stdv;
    //int32_t pos;   //todo : always -1 can be removed
    //int32_t state; //todo : always -1 can be removed
} event_t;

/* event table : adapted from scrappie */
typedef struct {
    size_t n;     //todo : int32_t not enough?
    size_t start; //todo : always 0?
    size_t end;   //todo : always equal to n?
    event_t* event;
} event_table;

/* k-mer model */
typedef struct {
    float level_mean;
    float level_stdv;

#ifdef CACHED_LOG
    float level_log_stdv;     //pre-calculated for efficiency
#endif

#ifdef LOAD_SD_MEANSSTDV
    //float sd_mean;
    //float sd_stdv;
    //float weight;
#endif
} model_t;

typedef struct {
    int32_t num_ref;
    char **ref_names;
    int32_t *ref_lengths;
    int32_t *ref_seq_lengths;

    float **forward;
    float **reverse;
} refsynth_t;

/* scaling parameters for the signal : taken from nanopolish */
typedef struct {
    // direct parameters that must be set
    float scale;
    float shift;
    //float drift; = 0 always?
    float var; // set later when calibrating
    //float scale_sd;
    //float var_sd;

#ifdef CACHED_LOG
    float log_var;    // derived parameters that are cached for efficiency
#endif
    //float scaled_var;
    //float log_scaled_var;
} scalings_t;

#endif
