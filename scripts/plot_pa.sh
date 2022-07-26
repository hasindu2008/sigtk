#!/bin/bash

#scripts/plot_pa.sh test/sequin_rna.blow5 00213403-4297-4f03-8412-3cc8b9cb845a

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <file.blow5> <read_id>"
    exit 1
fi

# Colour codes for printing
GREEN="\e[32m"
RED="\e[31m"
NORMAL="\033[0;39m"

[ -z ${SIGTK} ] && export SIGTK=sigtk
${SIGTK} --version &> /dev/null || { echo -e $RED"sigtk not found! Either put sigtk under path or set SIGTK variable, e.g.,export SIGTK=/path/to/sigtk"$NORMAL; exit 1;}

FILE=${1}
read_id=${2}

rm -f sigtk_${read_id}.pa.tmp

${SIGTK} pa -n "$FILE" ${read_id} |  cut -f 3 | tr ',' '\n' > sigtk_${read_id}.pa.tmp

if [[ "${SIGTK_PLOT_MTD}" == "matlab" ]]; then
    matlab.exe -nosplash -nodesktop -r "
    a=dlmread('sigtk_${read_id}.pa.tmp');
    plot(a); hold on; xlabel('sample index'), ylabel('current value (pA)'); legend('pA signal'); title('$read_id');
    "
elif [[ "${SIGTK_PLOT_MTD}" == "gnuplot" ]]; then
    gnuplot -e "
    set terminal pngcairo size 800,600 enhanced font 'Verdana,10';
    set output 'sigtk_${read_id}.pa.png';
    set xlabel 'sample index';
    set ylabel 'raw signal value';
    set title '$read_id';
    plot 'sigtk_${read_id}.pa.tmp' with lines title 'pA signal';
    "
    echo -e $GREEN"Plot written to sigtk_${read_id}.pa.png"$NORMAL
    rm -f sigtk_${read_id}.pa.tmp
else
    echo -e $RED"SIGTK_PLOT_MTD variable not set properly! Either set SIGTK_PLOT_MTD to matlab or gnuplot. e.g.,export SIGTK_PLOT_MTD=gnuplot"$NORMAL
    exit 1
fi



