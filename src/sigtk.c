/* @file sigtk.c
**
** @@
******************************************************************************/
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sigtk.h"
#include "misc.h"

#include <slow5/slow5.h>

#include <sys/wait.h>
#include <unistd.h>

enum sigtk_log_level_opt _log_level = LOG_VERB;

enum sigtk_log_level_opt get_log_level(){
    return _log_level;
}

void set_log_level(enum sigtk_log_level_opt level){
    _log_level = level;
}