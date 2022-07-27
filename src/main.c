/* @file main.c
**
******************************************************************************/
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <assert.h>
#include "sigtk.h"
#include "misc.h"
#include "error.h"

#ifdef HAVE_EXECINFO_H
    #include <execinfo.h>
#endif

extern enum sigtk_log_level_opt sigtk_log_level;

//make the segmentation faults a bit cool
void sig_handler(int sig) {
#ifdef HAVE_EXECINFO_H
    void* array[100];
    size_t size = backtrace(array, 100);
    ERROR("I regret to inform that a segmentation fault occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
    fprintf(stderr,
            "[%s::DEBUG]\033[1;35m Here is the backtrace in case it is of any "
            "use:\n",
            __func__);
    backtrace_symbols_fd(&array[2], size - 1, STDERR_FILENO);
    fprintf(stderr, "\033[0m\n");
#else
    ERROR("I regret to inform that a segmentation fault occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
#endif
    exit(EXIT_FAILURE);
}

int cmain(int argc, char* argv[], char *mode);
int srefmain(int argc, char* argv[]);

int print_usage(FILE *fp_help){

    fprintf(fp_help,"Usage: sigtk <command> [options]\n\n");
    fprintf(fp_help,"command:\n");
    fprintf(fp_help,"         sref      general synthetic signal for a reference\n");
    fprintf(fp_help,"         event     segment raw signal into events\n");
    fprintf(fp_help,"         stat      print statistics of the raw signal\n");
    // fprintf(fp_help,"         seg       segments such as adaptor and polyA\n");
    // fprintf(fp_help,"         jnn       print all segments found using james' neural network\n");
    fprintf(fp_help,"         pa        print raw signal in pico-amperes\n");
    if(fp_help==stderr){
        exit(EXIT_FAILURE);
    }
    else if(fp_help==stdout){
        exit(EXIT_SUCCESS);
    }
    else{
        assert(0);
    }


}

int main(int argc, char* argv[]){

    double realtime0 = realtime();
    signal(SIGSEGV, sig_handler);

    int ret=1;

    if(argc<2){
        return print_usage(stderr);
    }
    else if (strcmp(argv[1],"sref")==0){
        ret=srefmain(argc-1, argv+1);
    }
    else if (strcmp(argv[1],"event")==0 || strcmp(argv[1],"stat")==0 || strcmp(argv[1],"seg")==0 || strcmp(argv[1],"pa")==0 || strcmp(argv[1],"jnn")==0){
        ret=cmain(argc-1, argv+1, argv[1]);
    }
    else if(strcmp(argv[1],"--version")==0 || strcmp(argv[1],"-V")==0){
        fprintf(stdout,"sigtk %s\n",SIGTK_VERSION);
        exit(EXIT_SUCCESS);
    }
    else if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0){
        print_usage(stdout);
    }
    else{
        fprintf(stderr,"[sigtk] Unrecognised command %s\n",argv[1]);
        print_usage(stderr);
    }

    fprintf(stderr,"[%s] Version: %s\n", __func__,SIGTK_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }

    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec; Peak RAM: %.3f GB\n\n",
            __func__, realtime() - realtime0, cputime(),peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}
