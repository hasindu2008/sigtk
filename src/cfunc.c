/* @file  events_main.c
**
** @@
******************************************************************************/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "jnn.h"
#include "misc.h"
#include "sigtk.h"
#include "stat.h"

void print_events(char *rid, event_table et, slow5_rec_t *rec, opt_t opt){

    if(opt.compact){
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

void event_hdr(opt_t opt){
    if(opt.compact){
        printf("read_id\tlen_raw_signal\traw_start\traw_end\tnum_event\tevents\n");
    } else {
        printf("read_id\tevent_idx\traw_start\traw_end\tevent_mean\tevent_std\n");
    }
}


void event_func(slow5_rec_t *rec, opt_t opt){

    float *current_signal = signal_in_picoamps(rec);
    //trim(current_signal, rec->len_raw_signal);

    event_table et = getevents(rec->len_raw_signal, current_signal, opt.rna);

    print_events(rec->read_id,et, rec, opt);

    free(current_signal);
    free(et.event);
}

void pa_func(slow5_rec_t *rec, opt_t opt){

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

void jnn_func(slow5_rec_t *rec, opt_t opt){
    printf("%s\t",rec->read_id);

    uint64_t len_raw_signal = rec->len_raw_signal;
    printf("%ld\t",len_raw_signal);

    jnn_print(rec, opt.compact, opt.rna);
    printf("\n");

}
void jnn_hdr(){
    printf("read_id\tlen_raw_signal\tnum_seg\tseg\n");
}

void stat_hdr(){
    printf("read_id\tlen_raw_signal\traw_mean\tpa_mean\traw_std\tpa_std\traw_median\tpa_median\n");
}

void stat_func(slow5_rec_t *rec, opt_t opt){
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

void prefix_hdr(opt_t opt){
    printf("read_id\tlen_raw_signal\tadapt_start\tadapt_end\tpolya_start\tpolya_end");
    if(opt.p_stat){
        printf("\tadapt_mean\tadapt_std\tadapt_median\tpolya_mean\tpolya_std\tpolya_median");
    }
    printf("\n");
}

void prefix_func(slow5_rec_t *rec, opt_t opt){
    int64_t len_raw_signal = rec->len_raw_signal;
    printf("%s\t%ld\t",rec->read_id, len_raw_signal);
    jnn_pair_t p=find_adaptor(rec);
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


        jnn_pair_t polya;
        if(opt.rna){
            polya=find_polya(adapt_end,len_raw_signal-p.y, m_a+30+20,m_a+30-20);
        } else {
            polya.x = -1;
            polya.y = -1;
        }
        if(polya.y > 0){
            assert(polya.y + p.y < len_raw_signal);
            printf("%ld\t%ld",polya.x+p.y, polya.y+p.y);
        }else{
            printf(".\t.");
        }

        if(opt.p_stat){
            printf("\t%f\t%f\t%f\t",m_a, s_a, k_a);

            if(polya.y > 0){
                float m_poly = meanf(&current[polya.x+p.y],polya.y-polya.x);
                float s_poly = stdvf(&current[polya.x+p.y],polya.y-polya.x);
                float k_poly = medianf(&current[polya.x+p.y],polya.y-polya.x);
                printf("\t%f\t%f\t%f\t",m_poly, s_poly, k_poly);
            } else {
                printf("\t.\t.\t.");
            }
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
