#!/bin/bash
FASTQ_FILES=/home/anton/data/SRC/raw_data/*.fastq.gz
ADPTR_SHORT_5="GCTCTTCCGATCT"
ADPTR_SHORT_3="AGATCGGAAGAGC"
for fq in ${FASTQ_FILES}; do
	L=9
	fq_base=`basename ${fq}`
	while [ ${L} -le 13 ]; do
		cutadapt -g "${ADPTR_SHORT_5}" -a "${ADPTR_SHORT_3}" -n 2 -O ${L} --match-read-wildcards ${fq} -o /dev/null > clip_${L}_${fq_base%.fastq.gz}.stats
		grep "Trimmed reads" clip_${L}_${fq_base%.fastq.gz}.stats | sed -e 's/^..................................\(\)/\1/' | sed 's/..$//' >> clip_common_${fq_base%.fastq.gz}.stats
		rm clip_${L}_${fq_base%.fastq.gz}.stats
		let L=L+1
	done
	echo "I'm cut adapters from ${fq_base}"
done
