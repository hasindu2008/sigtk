/* @file  jnn.c
**
** @@
******************************************************************************/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "error.h"
#include "jnn.h"
#include "stat.h"

#define OUTLIER_MAX 1200
#define OUTLIER_MIN 0

static float *rolling_window(const float *x, int n, int w) {
    // int i = 0;

    assert(w<n);

    float *t = (float*)malloc(sizeof(float)*(n-w));
    MALLOC_CHK(t);

    // for(int i=0;i<n-w;i++){
    //     float tt = 0.0;
    //     for(int j=0;j<w;j++){
    //         tt += x[i+j];
    //     }
    //     t[i]=tt/w;
    // }


    float tt = 0.0;
    for(int i=0;i<w;i++){
        tt += x[i];
    }
    t[0]=tt/w;
    for(int i=1;i<n-w;i++){
        tt -= x[i-1];
        assert(i+w-1<n);
        tt += x[i+w-1];
        t[i]=tt/w;
    }


    // for(int i=0; i<n-w; i++){ //can remove later with a single var that inits to 0
    //     fprintf(stderr,"%f\t",t[i]);
    // }
    // fprintf(stderr,"\n");

    return t;
}


//adapted from
//https://github.com/Psy-Fer/deeplexicon/blob/master/scripts/dRNA_segmenter.py
static float *rm_outlier(const int16_t *x, int n) {

    float *t = (float*)malloc(sizeof(float)*(n));
    MALLOC_CHK(t);

    for(int i=0; i<n; i++){
        if(x[i]>OUTLIER_MAX){
            t[i] =OUTLIER_MAX;
        } else if (x[i]<OUTLIER_MIN){
            t[i] = OUTLIER_MIN;
        } else {
            t[i] = x[i];
        }
    }

    return t;
}

static float *rm_outlierf(const float *x, int n) {

    float *t = (float*)malloc(sizeof(float)*(n));
    MALLOC_CHK(t);

    for(int i=0; i<n; i++){
        if(x[i]>OUTLIER_MAX){
            t[i] =OUTLIER_MAX;
        } else if (x[i]<OUTLIER_MIN){
            t[i] = OUTLIER_MIN;
        } else {
            t[i] = x[i];
        }
    }

    return t;
}

//adapted from
//https://github.com/Psy-Fer/deeplexicon/blob/master/scripts/dRNA_segmenter.py
jnn_pair_t jnnv2(const int16_t *sig, int64_t nsample, jnnv2_param_t param){

    if(nsample > param.window){

        float *current = rm_outlier(sig,nsample);
        float *t = rolling_window(current,nsample,param.window);
        free(current);
        float mn = meanf(t,nsample-param.window);
        float std = stdvf(t,nsample-param.window);

        //int top = (MEAN_VAL + (STDV_VAL*0.5));
        float bot = mn - (std*param.std_scale);

        int8_t begin = 0;
        int seg_dist = param.seg_dist;
        int hi_thresh = param.hi_thresh;
        int lo_thresh = param.lo_thresh;
        int start=0;
        int end=0;

        int seg_c = SIGTK_SIZE;
        jnn_pair_t *segs = (jnn_pair_t *)malloc(sizeof(jnn_pair_t)*seg_c);
        int seg_i = 0;

        int count = -1;

        for (int j=0; j<nsample-param.window; j++){
            float i = t[j];
            count++;
            if (i < bot && !begin){
                start = count;
                begin = 1;
            }
            else if (i < bot){
                end = count;
            }
            else if (i > bot && begin){
                if (seg_i && start - segs[seg_i-1].y < seg_dist){
                    segs[seg_i-1].y = end;
                }
                else{
                    if(seg_i>=seg_c){
                        seg_c *= 2;
                        segs = (jnn_pair_t *)realloc(segs,sizeof(jnn_pair_t)*seg_c);
                    }
                    segs[seg_i].x = start;
                    segs[seg_i].y = end;
                    seg_i++;
                }
                start = 0;
                end = 0;
                begin = 0;
            }
        }
        jnn_pair_t p = {0, 0};
        for (int i=0; i<seg_i; i++){
            int b = segs[i].y;
            int a = segs[i].x;
            if (b - a > hi_thresh){
                continue;
            }
            if (b - a < lo_thresh){
                continue;
            }
            p.x = a+param.window/2-1;
            p.y = b+param.window/2-1;
            //fprintf(stderr,"FF %d\t%d\n",p.x, p.y);
            break;
        }

        free(t);
        free(segs);
        return p;
    }
    else{
        WARNING("%s","Not enough data to trim\n");
        jnn_pair_t p = {-1,-1};
        return p;
    }

}

jnn_pair_t find_adaptor(slow5_rec_t *rec){
    jnnv2_param_t param = JNNV2_RNA_ADAPTOR;
    return jnnv2(rec->raw_signal, rec->len_raw_signal, param);
}

jnn_pair_t *jnn_core(const float *sig, int64_t nsample, jnn_param_t param, int *n){

    float top;
    float bot;

    if(param.std_scale > 0){
        float mn = meanf(sig,nsample);
        float std = stdvf(sig,nsample);
        top = mn + (std*param.std_scale);
        bot = mn - (std*param.std_scale);
    } else {
        top = param.top;
        bot = param.bot;
    }

    int8_t prev = 0;    // previous string
    int err = 0;        // total error
    int prev_err = 0;  // consecutive error
    int c = 0;         // counter
    int w = param.corrector;   // window to increase total error thresh
    int seg_dist = param.seg_dist; // distance between 2 segs to be merged as one
    int start = 0;     // start pos
    int end = 0 ;      // end pos
    int window = param.window;
    int error = param.error;
    float stall_len = param.stall_len;

    // segments [(start, stop)]
    int seg_c = SIGTK_SIZE;
    jnn_pair_t * segs = (jnn_pair_t *)malloc(sizeof(jnn_pair_t)*seg_c);
    int seg_i = 0;

    for(int i=0; i<nsample; i++){
        float a = sig[i];
        if (a < top && a > bot){ // If datapoint is within range
            if (!prev){
                start = i;
                prev = 1;
            }
            c++; // increase counter
            w++; // increase window corrector count
            if (prev_err){
                prev_err = 0;
            }
            if (c >= window && c >= w &&  !(c % w)){ // if current window longer than detect limit, and corrector, and is divisible by corrector
                err--; // drop current error count by 1
            }
        }
        else{
            if (prev && err < error){
                c++;
                err++;
                prev_err++;
                if (c >= window && c >= w && !(c % w)){
                    err--;
                }
            }
            else if (prev && (c >= window || (!seg_i && c >= window * stall_len))){
                end = i - prev_err; // go back to where error stretch began for accurate cutting
                prev = 0;
                if (seg_i && start - segs[seg_i-1].y < seg_dist){ // if segs very close, merge them
                    segs[seg_i-1].y = end;
                }
                else{
                    if(seg_i>=seg_c){
                        seg_c *= 2;
                        segs = (jnn_pair_t *)realloc(segs,sizeof(jnn_pair_t)*seg_c);
                    }
                    segs[seg_i].x = start;
                    segs[seg_i].y = end;
                    seg_i++;
                }
                c = 0;
                err = 0;
                prev_err = 0;
            }
            else if (prev){
                prev = 0;
                c = 0;
                err = 0;
                prev_err = 0;
            }
        }
    }

    *n = seg_i;
    return segs;

}


//adapted from https://github.com/Psy-Fer/SquiggleKit/blob/a667e461b82f0ccc0a8103a8c1759515d6e34ac9/segmenter.py
jnn_pair_t *jnn_raw(const int16_t *raw, int64_t nsample, jnn_param_t param, int *n){

    jnn_pair_t *segs = NULL;
    *n = 0;
    if(nsample > 0){
        float *sig = rm_outlier(raw,nsample);
        segs = jnn_core(sig, nsample, param, n);
        free(sig);
    }

    return segs;
}

jnn_pair_t *jnn_pa(const float *raw, int64_t nsample, jnn_param_t param, int *n){

    jnn_pair_t *segs = NULL;
    *n = 0;
    if(nsample > 0){
        float *sig = rm_outlierf(raw,nsample);
        segs = jnn_core(sig, nsample, param, n);
        free(sig);
    }

    return segs;
}


jnn_pair_t jnn_print(slow5_rec_t *rec, int8_t fmt){

    int seg_i = 0;
    jnn_param_t param = JNNV1_PARAM;
    jnn_pair_t *segs = jnn_raw(rec->raw_signal,rec->len_raw_signal,param,&seg_i);
    jnn_pair_t p = {-1,-1};

    if(segs){
        printf("%d\t",seg_i);

        if(fmt){
            uint64_t ci = 0;
            uint64_t mi = 0;
            for(int i=0; i<seg_i; i++){
                ci += (mi = segs[i].x - ci);
                if(mi) printf("%dH",(int)mi);
                ci += (mi = segs[i].y - ci);
                if(mi) printf("%d,",(int)mi);
            }
        } else {
            for(int i=0; i<seg_i; i++){
                printf("%ld,%ld;",segs[i].x,segs[i].y);
            }
        }
        if(seg_i){
            p = segs[seg_i-1];
        } else {
            printf(".");
        }
        free(segs);
    }

    return p;

}

jnn_pair_t find_polya(const float *raw, int64_t nsample, float top, float bot){
    jnn_pair_t p = {-1,-1};
    int seg_i = 0;
    jnn_param_t param = JNNV1_POLYA;
    param.top = top;
    param.bot = bot;
    jnn_pair_t *segs = jnn_pa(raw,nsample,param,&seg_i);

    if(segs){
        if(seg_i > 0){
            p = segs[0];
        }
        free(segs);
    }


    return p;
}
