#!/bin/bash
DIR=/home/anton/data/debug/temp
fq=$1
prefix=$2
cd $DIR
mkdir $prefix
count=5
FASTX_REVCOM=fastx_reverse_complement
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #CTCGGCCGCGACC
cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $DIR/$prefix/inner0-gatcs.fastq $DIR/$fq -o $DIR/$prefix/output0-gatc.fastq >> $DIR/$prefix/clip0_gatc.stats
while [ $count -le 13 ]; do
count_for_cut=$((13-$((${count}-4))))
count_for_name=$((${count}-4))
adptr5=`echo $ADPTR_SHORT_5 | sed -e 's/^.\{'$count_for_cut'\}//'`
adptr3=`echo $ADPTR_SHORT_3 | sed -e 's/.\{'$count_for_cut'\}$//'`
cutadapt -g "^${adptr5}GATC" -a "GATC${adptr3}$" -O $count -e 0.01 --no-trim --untrimmed-output $DIR/$prefix/inner${count_for_name}-gatcs.fastq $DIR/$fq -o $DIR/$prefix/output${count_for_name}-gatc.fastq >> $DIR/$prefix/clip${count_for_name}-gatc.stats
count=$(($count+1))
done
