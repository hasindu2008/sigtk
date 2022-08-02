/* @file misc.h
**
** jnn
** @@
******************************************************************************/

#ifndef JNN_H
#define JNN_H

#include <stdint.h>
#include <slow5/slow5.h>

typedef struct {
    int64_t x;
    int64_t y;
} jnn_pair_t;

typedef struct {
    float std_scale;
    int corrector; //corrector, window to increase total error thresh
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    float stall_len;
    int error;
    float top; //will be only used if std_scale is -1
    float bot; //will be only used if std_scale is -1
} jnn_param_t;


#define JNNV1_PARAM { \
    .std_scale = 0.75, \
    .corrector = 50, \
    .seg_dist = 50, \
    .window = 1000, \
    .stall_len = 1.0, \
    .error = 5, \
    .top = 0, \
    .bot = 0, \
} \

//relative
#define JNNV1_POLYA { \
    .std_scale = -1, \
    .corrector = 50, \
    .seg_dist = 200, \
    .window = 250, \
    .stall_len = 1.0, \
    .error = 30, \
    .top = 0, \
    .bot = 0, \
} \


typedef struct {
    float std_scale;
    int seg_dist; // distance between 2 segs to be merged as one
    int window;
    float stall_len;
    int hi_thresh;
    int lo_thresh;
} jnnv2_param_t;

//dRNA segmenter
#define JNNV2_RNA_ADAPTOR { \
    .std_scale = 0.5, \
    .seg_dist = 1500, \
    .window = 2000, \
    .hi_thresh = 200000, \
    .lo_thresh = 2000, \
} \

#define SIGTK_SIZE 1000

jnn_pair_t *jnn_raw(const int16_t *raw, int64_t nsample, jnn_param_t param, int *n);
jnn_pair_t *jnn_pa(const float *raw, int64_t nsample, jnn_param_t param, int *n);
jnn_pair_t find_polya(const float *raw, int64_t nsample, float top, float bot);
jnn_pair_t find_adaptor(slow5_rec_t *rec);
jnn_pair_t jnn_print(slow5_rec_t *rec, int8_t fmt);
jnn_pair_t jnnv2(const int16_t *sig, int64_t nsample, jnnv2_param_t param);
#endif
