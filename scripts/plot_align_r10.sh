#!/bin/bash

#This is a dirty script to plot the alignment of a signal to a reference using matlab
#mainly for my own use

#scripts/plot_align.sh test/batch0.blow5 1d72c1af-0e12-4326-9310-407a15fa7b2b test/batch0.fastq

set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 <file.blow5> <read_id> <reads.fastq>"
    exit 1
fi

# Colour codes for printing
GREEN="\e[32m"
RED="\e[31m"
NORMAL="\033[0;39m"

FILE=${1}
read_id=${2}
fastq=${3}

[ -z ${SIGTK} ] && export SIGTK=sigtk
${SIGTK} --version &> /dev/null || { echo -e $RED"sigtk not found! Either put sigtk under path or set SIGTK variable, e.g.,export SIGTK=/path/to/sigtk"$NORMAL; exit 1;}
[ -z ${SLOW5TOOLS} ] && export SLOW5TOOLS=slow5tools
${SLOW5TOOLS} --version &> /dev/null || { echo -e $RED"slow5tools not found! Either put slow5tools under path or set SLOW5TOOLS variable, e.g.,export SLOW5TOOLS=/path/to/slow5tools"$NORMAL; exit 1;}
[ -z ${F5C} ] && export F5C=f5c
${F5C} --version &> /dev/null || { echo -e $RED"f5c not found! Either put f5c under path or set F5C variable, e.g.,export F5C=/path/to/f5c"$NORMAL; exit 1;}
samtools --version &> /dev/null || { echo -e $RED"samtools not found in path!"$NORMAL; exit 1;}

rm -f sigtk_${read_id}.tmp sigtk_${read_id}.events.tmp sigtk_${read_id}.fasta sigtk_${read_id}.blow5 sigtk_${read_id}.resquigged.tmp sigtk_${read_id}_bases.tmp

KMER_SIZE=9
${SLOW5TOOLS} get ${FILE} ${read_id} -o sigtk_${read_id}.blow5 || { echo -e $RED"Error: slow5tools get failed"$NORMAL; exit 1;}
${SLOW5TOOLS} view sigtk_${read_id}.blow5 | grep "^[@#]" | grep "experiment_type" | grep "rna" && KMER_SIZE=5;
echo "KMER_SIZE: ${KMER_SIZE}"
${SLOW5TOOLS} get --to slow5 ${FILE} ${read_id} | grep -v '^[#@]' | awk '{print $8}' > sigtk_${read_id}.tmp || { echo -e $RED"Error: slow5tools get failed"$NORMAL; exit 1;}
${SIGTK} event ${FILE} ${read_id} -n -c | awk '{print $3"\t"$6}' | tr ',' '\t' > sigtk_${read_id}.events.tmp || { echo -e $RED"Error: sigtk event failed"$NORMAL; exit 1;}
samtools faidx ${fastq} ${read_id} > sigtk_${read_id}.fasta || { echo -e $RED"Error: samtools faidx failed"$NORMAL; exit 1;}
${F5C} resquiggle sigtk_${read_id}.fasta sigtk_${read_id}.blow5 --kmer-model "$(dirname $(which f5c))/test/r10-models/r10.4.1_400bps.nucleotide.9mer.template.model" | tail -n+2 | awk '{print $3"\t"$4}' | sed 's/\./-1/g' > sigtk_${read_id}.resquigged.tmp || { echo -e $RED"Error: f5c resquiggle failed"$NORMAL; exit 1;}
grep -v '^>' sigtk_${read_id}.fasta | tr -d '\n' > sigtk_${read_id}_bases.tmp || { echo -e $RED"Error: Read not found in fastq"$NORMAL; exit 1;}
rm -f sigtk_${read_id}.fasta sigtk_${read_id}.blow5

if [[ "${SIGTK_PLOT_MTD}" == "matlab" ]]; then

    matlab.exe -nosplash -nodesktop -minimize -r "
    a=dlmread('sigtk_${read_id}.tmp'); b=dlmread('sigtk_${read_id}.events.tmp');
    %startidx=b(:,1)+1; endidx=b(:,2);
    %avg=zeros(length(a),1);
    %for j=1:length(startidx)
    %    avg(startidx(j):endidx(j))=mean(a(startidx(j):endidx(j)));
    %end

    startidx=b(1)+1; % b(1) is raw start, +1 for 1-based indexing
    n_events=length(b)-1
    avg=zeros(length(a),1);
    for j=1:n_events
        endidx=startidx+b(j+1)-1; % make closed interval
        avg(startidx:endidx)=mean(a(startidx:endidx));
        startidx=endidx+1;
    end

    cmap = colororder();
    rsq = dlmread('sigtk_${read_id}.resquigged.tmp');
    base = readlines('sigtk_${read_id}_bases.tmp'); c=cellstr(base);
    startidx2=rsq(:,1)+1;
    endidx2=rsq(:,2);
    avg2=zeros(length(a)-${KMER_SIZE}+1,1);
    x=zeros(strlength(base)-${KMER_SIZE}+1,1);
    y=zeros(strlength(base)-${KMER_SIZE}+1,1);
    if (strlength(base)-${KMER_SIZE}+1 ~= length(startidx2))
        disp('Error: lengths do not match');
    end
    prevx=0; prevy=0;
    for j=1:length(startidx2)
        if(startidx2(j)>0 && endidx2(j)>0)
            avg2(startidx2(j):endidx2(j))=mean(a(startidx2(j):endidx2(j)));
            x(j)=(startidx2(j)+endidx2(j))/2;
            y(j)=min(a(startidx2(j)),a(endidx2(j)));
            prevx=x(j); prevy=y(j);
        else
            x(j)=prevx;
            y(j)=prevy+100;
        end
    end

    plot(a,'Color',cmap(1,:)); hold on; plot(avg,'Color',cmap(3,:), 'LineStyle', '--'); plot(avg2,'Color',cmap(2,:)), xlabel('sample index'), ylabel('raw signal value');  legend('raw signal','events','align');

    %kmers
    base=convertStringsToChars(base)';

    for j=1:length(base)-${KMER_SIZE}+1
        if(${KMER_SIZE}==9)
            kmers{j}=strcat(char(base(j)),char(base(j+1)),char(base(j+2)),char(base(j+3)),char(base(j+4)),char(base(j+5)),char(base(j+6)),char(base(j+7)),char(base(j+8)));
        else
            kmers{j}=strcat(char(base(j)),char(base(j+1)),char(base(j+2)),char(base(j+3)),char(base(j+4)));
        end
    end

    if(${KMER_SIZE}==5)
        kmers=flip(kmers);
    end
    h=text(x,y,cellstr(kmers),'FontSize',7);
    set(h,'Rotation',90);

    "
else
    echo -e $RED"SIGTK_PLOT_MTD variable not set properly! set SIGTK_PLOT_MTD to matlab. e.g.,export SIGTK_PLOT_MTD=matlab"$NORMAL
    exit 1
fi

