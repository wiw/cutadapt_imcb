#!/bin/bash
FASTQ_FILES=/home/anton/data/SRC/2047_3_1.3.1_Dam_Scrb__1.3.1_Dam_Scrb__ACAGTGA_L003_R1_001.fastq.gz
FASTX_REVCOM=fastx_reverse_complement
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
for fq in ${FASTQ_FILES}; do
	L=13
	fq_base=`basename ${fq}`
	while [ ${L} -le 13 ]; do
		cutadapt -g "${ADPTR_SHORT_5}" -a "${ADPTR_SHORT_3}" -n 2 -O ${L} --match-read-wildcards ${fq} -o output > clip_${L}_${fq_base%.fastq.gz}.stats | cutadapt -g "${ADPTR_SHORT_5}" -a "${ADPTR_SHORT_3}" -n 2 -O ${L} --match-read-wildcards ${fq} -o /dev/null > clip_${L}_${fq_base%.fastq.gz}.stats
		grep "Trimmed reads" clip_${L}_${fq_base%.fastq.gz}.stats | sed -e 's/^..................................\(\)/\1/' | sed 's/..$//' >> clip_common_${fq_base%.fastq.gz}.stats
		rm clip_${L}_${fq_base%.fastq.gz}.stats
		let L=L+1
	done
done
