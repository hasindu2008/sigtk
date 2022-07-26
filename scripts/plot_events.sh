#!/bin/bash

#scripts/plot_events.sh test/sequin_rna.blow5 00213403-4297-4f03-8412-3cc8b9cb845a

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

rm -f sigtk_${read_id}.tmp sigtk_${read_id}.events.tmp

slow5tools get --to slow5 ${FILE} ${read_id} | grep -v '^[#@]' | awk '{print $8}' > sigtk_${read_id}.tmp
sigtk event ${FILE} ${read_id} -n | awk '{print $3"\t"$4"\t"$5}' > sigtk_${read_id}.events.tmp

if [[ "${SIGTK_PLOT_MTD}" == "matlab" ]]; then

    matlab.exe -nosplash -nodesktop -minimize -r "
    a=dlmread('sigtk_${read_id}.tmp'); b=dlmread('sigtk_${read_id}.events.tmp');
    startidx=b(:,1)+1; endidx=b(:,2);
    avg=zeros(length(a),1);
    for j=1:length(startidx)
        avg(startidx(j):endidx(j))=mean(a(startidx(j):endidx(j)));
    end
    plot(a); hold on; plot(avg); xlabel('sample index'), ylabel('raw signal value');  legend('raw signal','events');
    figure; plot(b(:,3));
    "
else
    echo -e $RED"SIGTK_PLOT_MTD variable not set properly! set SIGTK_PLOT_MTD to matlab. e.g.,export SIGTK_PLOT_MTD=matlab"$NORMAL
    exit 1
fi

