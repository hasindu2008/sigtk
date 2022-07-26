/* @file  events_main.c
**
** @@
******************************************************************************/
#define _XOPEN_SOURCE 700
#include "sigtk.h"
#include "misc.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

extern int compact;

float ks_ksmall_float(size_t, float*, size_t);
int16_t ks_ksmall_int16_t(size_t, int16_t*, size_t);

void print_events(char *rid, event_table et, slow5_rec_t *rec){

    if(compact){
        printf("%s\t%ld\t", rid, rec->len_raw_signal);
        if(et.n){
            assert(et.event[et.n-1].start < rec->len_raw_signal);
            assert(et.event[et.n-1].start+(int)et.event[et.n-1].length <= rec->len_raw_signal);
            printf("%ld\t%ld\t", et.event[0].start, et.event[et.n-1].start+(int)et.event[et.n-1].length);
            printf("%ld\t", et.n);
            uint32_t j = 0;
            uint64_t ci = et.event[0].start; //current index
            uint64_t mi = 0;
            for (j = 0; j < et.n; j++) {
                ci += (mi =  et.event[j].start - ci);
                //if(mi) printf("%d;",(int)mi);
                assert(mi==0); //do not support jumps at the moment
                ci += (mi = et.event[j].length);
                if(mi) {
                    if (j<et.n-1) printf("%d,",(int)mi); else printf("%d",(int)mi);
                }
            }
        }
        else{
            printf(".\t.\t.\t.");
        }


        //printf("\t");

        // for (j = 0; j < et.n; j++) {
        //     printf("%ld,%d,%f,%f;", et.event[j].start,
        //             (int)et.event[j].length, et.event[j].mean,
        //             et.event[j].stdv);
        // }
        printf("\n");
    } else {
        uint32_t j = 0;
        for (j = 0; j < et.n; j++) {
            printf("%s\t%d\t%ld\t%ld\t%f\t%f\n", rid, j, et.event[j].start,
                    et.event[j].start+(int)et.event[j].length, et.event[j].mean,
                    et.event[j].stdv);
        }
        printf("\n");
    }

}

void event_hdr(){
    if(compact){
        printf("read_id\tlen_raw_signal\traw_start\traw_end\tnum_event\tevents\n");
    } else {
        printf("read_id\tevent_idx\traw_start\traw_end\tevent_mean\tevent_std\n");
    }
}


void event_func(slow5_rec_t *rec, int8_t rna){

    float *current_signal = signal_in_picoamps(rec);
    //trim(current_signal, rec->len_raw_signal);

    event_table et = getevents(rec->len_raw_signal, current_signal, rna);


    print_events(rec->read_id,et, rec);

    free(current_signal);
    free(et.event);
}

void pa_func(slow5_rec_t *rec, int8_t rna){

    float *current_signal = signal_in_picoamps(rec);

    uint64_t len_raw_signal = rec->len_raw_signal;
    printf("%s\t%ld\t",rec->read_id,len_raw_signal);

    for(uint64_t i=0;i<len_raw_signal;i++){
        if (i==len_raw_signal-1) {
            printf("%f",current_signal[i]);
        } else {
            printf("%f,",current_signal[i]);
        }
    }
    printf("\n");

    free(current_signal);
}

void pa_hdr(){
    printf("read_id\tlen_raw_signal\tpa\n");
}

float meanf(float *x, int n) {
    float sum = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}

float meani16(int16_t *x, int n) {
    float sum = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}


float stdvf(float *x, int n) {
    float sum = 0;
    int i = 0;
    float m = meanf(x, n); //can reuse already calculated mean
    for (i = 0; i < n; i++) {
        sum += (x[i] - m) * (x[i] - m);
    }
    return sqrtf(sum / n);
}

float stdvi16(int16_t *x, int n) {
    float sum = 0;
    int i = 0;
    float m = meani16(x, n); //can reuse already calculated mean
    for (i = 0; i < n; i++) {
        sum += (x[i] - m) * (x[i] - m);
    }
    return sqrtf(sum / n);
}

float medianf(float* x, int n) {

    float *copy = (float *)malloc(n * sizeof(float));
    memcpy(copy, x, n * sizeof(float));
    float m = ks_ksmall_float(n, copy, n / 2);
    free(copy);
    return m;

}

int16_t mediani16(int16_t* x, int n) {

    int16_t *copy = (int16_t *)malloc(n * sizeof(int16_t));
    memcpy(copy, x, n * sizeof(int16_t));
    int16_t m = ks_ksmall_int16_t(n, copy, n / 2);
    free(copy);
    return m;
}


#define OUTLIER_MAX 1200
#define OUTLIER_MIN 0

float *rolling_window(float *x, int n, int w) {
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
float *rm_outlier(int16_t *x, int n) {

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

float *rm_outlierf(float *x, int n) {

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
pair_t find_adaptor(slow5_rec_t *rec){

    int64_t nsample = rec->len_raw_signal;

    if(nsample > SIGTK_WINDOW_SIZE){


        float *current = rm_outlier(rec->raw_signal,rec->len_raw_signal);
        float *t = rolling_window(current,nsample,SIGTK_WINDOW_SIZE);
        free(current);
        float mn = meanf(t,nsample-SIGTK_WINDOW_SIZE);
        float std = stdvf(t,nsample-SIGTK_WINDOW_SIZE);

        //int top = (MEAN_VAL + (STDV_VAL*0.5));
        float bot = mn - (std*0.5);

        int8_t begin = 0;
        int seg_dist = 1500;
        int hi_thresh = 200000;
        int lo_thresh = 2000;
        int start=0;
        int end=0;

        int seg_c = SIGTK_SIZE;
        pair_t *segs = (pair_t *)malloc(sizeof(pair_t)*seg_c);
        int seg_i = 0;

        int count = -1;

        for (int j=0; j<nsample-SIGTK_WINDOW_SIZE; j++){
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
                        segs = (pair_t *)realloc(segs,sizeof(pair_t)*seg_c);
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
        pair_t p = {0, 0};
        for (int i=0; i<seg_i; i++){
            int b = segs[i].y;
            int a = segs[i].x;
            if (b - a > hi_thresh){
                continue;
            }
            if (b - a < lo_thresh){
                continue;
            }
            p.x = a+SIGTK_WINDOW_SIZE/2-1;
            p.y = b+SIGTK_WINDOW_SIZE/2-1;
            //fprintf(stderr,"FF %d\t%d\n",p.x, p.y);
            break;
        }

        free(t);
        free(segs);
        return p;
    }
    else{
        WARNING("%s","Not enough data to trim\n");
        pair_t p = {-1,-1};
        return p;
    }

}

//adapted from https://github.com/Psy-Fer/SquiggleKit/blob/a667e461b82f0ccc0a8103a8c1759515d6e34ac9/segmenter.py
pair_t jnn(slow5_rec_t *rec){

    int64_t nsample = rec->len_raw_signal;
    pair_t  p = {-1,-1};

    if(nsample > 0){

        float *sig = rm_outlier(rec->raw_signal,rec->len_raw_signal);

        float mn = meanf(sig,nsample);
        float std = stdvf(sig,nsample);

        float top = mn + (std*0.75);
        float bot = mn - (std*0.75);

        int8_t prev = 0;    // previous string
        int err = 0;        // total error
        int prev_err = 0;  // consecutive error
        int c = 0;         // counter
        int w = 50;        // window to increase total error thresh
        int seg_dist = 50; // distance between 2 segs to be merged as one
        int start = 0;     // start pos
        int end = 0 ;      // end pos
        int window = 1000;
        int error = 5;
        float stall_len = 1.0;

        // segments [(start, stop)]
        pair_t segs[SIGTK_SIZE];
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

        free(sig);

        printf("%d\t",seg_i);

        uint64_t ci = 0;
        uint64_t mi = 0;
        for(int i=0; i<seg_i; i++){
            ci += (mi = segs[i].x - ci);
            if(mi) printf("%dj",(int)mi);
            ci += (mi = segs[i].y - ci);
            if(mi) printf("%dm",(int)mi);
        }
        printf("\t");
        for(int i=0; i<seg_i; i++){
            printf("%ld,%ld;",segs[i].x,segs[i].y);
        }

        if(seg_i){
            p = segs[seg_i-1];
        }

    }

    return p;
}

void jnn_func(slow5_rec_t *rec, int8_t rna){
    printf("%s\t",rec->read_id);

    uint64_t len_raw_signal = rec->len_raw_signal;
    printf("%ld\t",len_raw_signal);

    jnn(rec);
    printf("\n");

}
void jnn_hdr(){
    printf("read_id\tlen_raw_signal\tnum_seg\tLSAR\tseg_st0,seg_end0;seg_st1,seg_end1;....\n");
}

pair_t find_polya(float *raw, int64_t nsample, float top, float bot){
    pair_t p = {-1,-1};
    if(nsample > 0){

        float *sig = rm_outlierf(raw,nsample);

        // float mn = mean(sig,nsample);
        // float std = stdv(sig,nsample);

        // float top = mn + (std*0.75);
        // float bot = mn - (std*0.75);

        int8_t prev = 0;    // previous string
        int err = 0;        // total error
        int prev_err = 0;  // consecutive error
        int c = 0;         // counter
        int w = 50;        // window to increase total error thresh
        int seg_dist = 200; // distance between 2 segs to be merged as one
        int start = 0;     // start pos
        int end = 0 ;      // end pos
        int window = 250;
        int error = 30;
        float stall_len = 1.0;

        // segments [(start, stop)]
        int seg_c = SIGTK_SIZE;
        pair_t *segs = (pair_t *)malloc(sizeof(pair_t)*seg_c);
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
                            segs = (pair_t *)realloc(segs,sizeof(pair_t)*seg_c);
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

        free(sig);

        // fprintf(stderr,"%d\t",seg_i);
        // for(int i=0; i<seg_i; i++){
        //     if (i==seg_i-1){
        //         fprintf(stderr,"%ld,%ld",segs[i].x,segs[i].y);
        //     }else{
        //         fprintf(stderr,"%ld,%ld\t",segs[i].x,segs[i].y);
        //     }
        // }
        // fprintf(stderr,"\n");

        if(seg_i){
            p = segs[0];
        }

        free(segs);

    }
    return p;
}

void stat_hdr(){
    printf("read_id\tlen_raw_signal\traw_mean\tpa_mean\traw_std\tpa_std\traw_median\tpa_median\n");
}

void stat_func(slow5_rec_t *rec, int8_t rna){
    printf("%s\t",rec->read_id);

    uint64_t len_raw_signal = rec->len_raw_signal;
    printf("%ld\t",len_raw_signal);

    float m1 = meani16(rec->raw_signal,len_raw_signal);
    float s1 = stdvi16(rec->raw_signal,len_raw_signal);
    int16_t k1 = mediani16(rec->raw_signal,len_raw_signal);

    float *current = signal_in_picoamps(rec);
    float m2 = meanf(current,len_raw_signal);
    float s2 = stdvf(current,len_raw_signal);
    float k2 = medianf(current,len_raw_signal);

    printf("%f\t%f\t%f\t%f\t%d\t%f\t",m1,m2,s1,s2,k1,k2);

    // float *t = rolling_window(current,len_raw_signal,SIGTK_WINDOW_SIZE);
    // int t_size = len_raw_signal-SIGTK_WINDOW_SIZE;

    // float m_t = mean(t,t_size);
    // float s_t = stdv(t,t_size);
    //

    // for(uint64_t i=0;i<len_raw_signal;i++){
    //     double pA = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
    //     printf("%f ",pA);
    // }


    // free(t);
    free(current);
    printf("\n");
}

void seg_hdr(){
    printf("read_id\tlen_raw_signal\tadapt_start\tadapt_end\tpolya_start\tpolya_end\tadapt_mean\tadapt_std\tadapt_median\tpolya_mean\tpolya_std\tpolya_median\n");
}

void seg_func(slow5_rec_t *rec, int8_t rna){
    int64_t len_raw_signal = rec->len_raw_signal;
    printf("%s\t%ld\t",rec->read_id, len_raw_signal);
    pair_t p=find_adaptor(rec);
    if(p.y > 0){
        assert(p.y<len_raw_signal);
        printf("%ld\t%ld\t",p.x, p.y);

        float *current = signal_in_picoamps(rec);
        float m_a = meanf(&current[p.x],p.y-p.x);
        float s_a = stdvf(&current[p.x],p.y-p.x);
        float k_a = medianf(&current[p.x],p.y-p.x);

        assert(p.y > 0);
        assert(p.y < len_raw_signal);


        float *adapt_end = &current[p.y];


        pair_t polya;
        if(rna){
            polya=find_polya(adapt_end,len_raw_signal-p.y, m_a+30+20,m_a+30-20);
        } else {
            polya.x = -1;
            polya.y = -1;
        }
        if(polya.y > 0){
            assert(polya.y + p.y < len_raw_signal);
            printf("%ld\t%ld\t",polya.x+p.y, polya.y+p.y);
        }else{
            printf(".\t.\t");
        }

        printf("%f\t%f\t%f\t",m_a, s_a, k_a);

        if(polya.y > 0){
            float m_poly = meanf(&current[polya.x+p.y],polya.y-polya.x);
            float s_poly = stdvf(&current[polya.x+p.y],polya.y-polya.x);
            float k_poly = medianf(&current[polya.x+p.y],polya.y-polya.x);
            printf("%f\t%f\t%f\t",m_poly, s_poly, k_poly);
        } else {
            printf(".\t.\t.");
        }

        free(current);

    } else{
        printf(".\t.\t.\t.");
    }




    // polya.x += p.y;
    // polya.y += p.y;
    // assert(polya.y < len_raw_signal);

    printf("\n");



    //trim_polya(rec);
}
