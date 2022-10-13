/* @file main.c
**
******************************************************************************/

#include <stdio.h>
#include <slow5/slow5.h>


int repmain(int argc, char* argv[]){

    if(argc < 2){
        fprintf(stderr,"Usage: %s <file>\n", argv[0]);
        return 1;
    }

    printf("sample\tcount_1\tcount_2\tcount_3\tcount_4\tcount_5\tcount_6\tcount_7\tcount_8\tcount_9+\n");

    for(int sample = 300; sample < 700; sample++){


        slow5_file_t *sp = slow5_open(argv[1],"r");
        if(sp==NULL){
        fprintf(stderr,"Error in opening file\n");
        return 1;
        }
        slow5_rec_t *rec = NULL;
        int ret=0;

        int64_t count_gg[10] = {0};

        while((ret = slow5_get_next(&rec,sp)) >= 0){

            uint64_t len_raw_signal = rec->len_raw_signal;
            int64_t count_g[10] = {0};

            int mark = 0;
            int count = 0;
            int c2 = 0;
            int c3 = 0;
            for(uint64_t i=0;i<len_raw_signal;i++){
                if (rec->raw_signal[i] == 455){
                    mark = 1;
                    count++;
                    c2++;
                } else {
                    if(mark){
                        if(count <10 ){
                            count_g[count]++;
                        } else {
                            count_g[9]++;
                        }
                    }
                    mark = 0;
                    count =0;
                }
            }
            if(mark){
                if(count <10 ){
                    count_g[count]++;
                } else {
                    count_g[9]++;
                }
            }

            for(int i=1;i<10;i++){
                count_gg[i] += count_g[i];
                c3 += count_g[i]*i;
            }
            if(count_g[9]==0 && c2!=c3){
                fprintf(stderr,"ASSERTION FAIL: %d\t%d\n",c2,c3);
                exit(EXIT_FAILURE);
            } else {
                if(c2<c3){
                    fprintf(stderr,"ASSERTION2 FAIL: %d\t%d\n",c2,c3);
                    exit(EXIT_FAILURE);

                }
            }
        }

        if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
            fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
            return 1;
        }

        slow5_rec_free(rec);

        slow5_close(sp);

        printf("%d\t",sample);
        for(int i=1;i<10;i++){
            if (i==9) {
                printf("%ld",count_gg[i]);
            } else {
                printf("%ld\t",count_gg[i]);
            }
        }
        printf("\n");

    }

    return 0;
}
