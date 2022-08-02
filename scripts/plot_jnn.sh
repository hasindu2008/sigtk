#!/bin/bash

#scripts/plot_jnn.sh test/sequin_rna.blow5 00213403-4297-4f03-8412-3cc8b9cb845a

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <file.blow5> <read_id>"
    exit 1
fi

# Colour codes for printing
GREEN="\e[32m"
RED="\e[31m"
NORMAL="\033[0;39m"

FILE=${1}
read_id=${2}

[ -z ${SIGTK} ] && export SIGTK=sigtk
${SIGTK} --version &> /dev/null || { echo -e $RED"sigtk not found! Either put sigtk under path or set SIGTK variable, e.g.,export SIGTK=/path/to/sigtk"$NORMAL; exit 1;}
[ -z ${SLOW5TOOLS} ] && export SLOW5TOOLS=slow5tools
${SLOW5TOOLS} --version &> /dev/null || { echo -e $RED"slow5tools not found! Either put slow5tools under path or set SLOW5TOOLS variable, e.g.,export SLOW5TOOLS=/path/to/slow5tools"$NORMAL; exit 1;}

rm -f sigtk_${read_id}.tmp sigtk_${read_id}.jnn.tmp

slow5tools get --to slow5 ${FILE} ${read_id} | grep -v '^[#@]' | awk '{print $8}' > sigtk_${read_id}.tmp || { echo -e $RED"Error: failed to get read_id ${read_id} from ${FILE}"$NORMAL; exit 1;}
sigtk jnn ${FILE} ${read_id} -n | cut -f 3,4 | tr ',;' '\t'  > sigtk_${read_id}.jnn.tmp || { echo -e $RED"Error: failed to apply jnn for read_id ${read_id} from ${FILE}"$NORMAL; exit 1;}
COUNT=$(cut -f1  sigtk_${read_id}.jnn.tmp);
[ "$COUNT" == "0" ] && { echo -e $RED"Error: no segments found for read_id ${read_id} from ${FILE}"$NORMAL; exit 1;}

if [[ "${SIGTK_PLOT_MTD}" == "matlab" ]]; then
    matlab.exe -nosplash -nodesktop -minimize -r "
    a=dlmread('sigtk_${read_id}.tmp'); b=dlmread('sigtk_${read_id}.jnn.tmp');
    c=b(:,1);
    plot(a); hold on;
    if(c>0)
        x=b(:,2:end);
        y=zeros(size(x))+1200;
        stem(x,y);
    end
    xlabel('sample index'), ylabel('raw signal value');  legend('raw signal','jnn'); title('$read_id');
    "
else
    echo -e $RED"SIGTK_PLOT_MTD variable not set properly! set SIGTK_PLOT_MTD to matlab. e.g.,export SIGTK_PLOT_MTD=matlab"$NORMAL
    exit 1
fi



