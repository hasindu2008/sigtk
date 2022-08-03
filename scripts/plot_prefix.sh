#!/bin/bash

#scripts/plot_prefix.sh test/sequin_rna.blow5 00213403-4297-4f03-8412-3cc8b9cb845a

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

rm -f sigtk_${read_id}.tmp sigtk_${read_id}.prefix.tmp

${SLOW5TOOLS} get --to slow5 ${FILE} ${read_id} | grep -v '^[#@]' | awk '{print $8}' > sigtk_${read_id}.tmp || { echo -e $RED"Error: failed to get read_id ${read_id} from ${FILE}"$NORMAL; exit 1;}
${SIGTK} prefix ${FILE} ${read_id} -n | cut -f 3,4,5,6 | sed 's/\./-1/g'  > sigtk_${read_id}.prefix.tmp || { echo -e $RED"Error: failed to get segments for read_id ${read_id} from ${FILE}"$NORMAL; exit 1;}

if [[ "${SIGTK_PLOT_MTD}" == "matlab" ]]; then

    matlab.exe -nosplash -nodesktop -minimize -r "
    a=dlmread('sigtk_${read_id}.tmp'); x=dlmread('sigtk_${read_id}.prefix.tmp');
    x1=x(:,[1:2]);
    x2=x(:,[3:4]);
    y=[1200,1200];
    plot(a); hold on;
    if(x1(1,1)>=0 && x1(1,2)>=0)
        stem(x1,y);
    end
    if(x2(1,1)>=0 && x2(1,2)>=0)
        stem(x2,y);
    end
    xlabel('sample index'), ylabel('raw signal value');
    legend('raw signal', 'adaptor','polya'); title('$read_id');
    "
else
    echo -e $RED"SIGTK_PLOT_MTD variable not set properly! set SIGTK_PLOT_MTD to matlab. e.g.,export SIGTK_PLOT_MTD=matlab"$NORMAL
    exit 1
fi




