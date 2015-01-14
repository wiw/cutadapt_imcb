#!/bin/bash
dir=/home/anton/Doc/Documents/SequenceData/DamID_IBprot/export_not_combine/141023
bzip_files=$dir/*.fastq.bz2
cd $dir
for fq in $bzip_files; do
	fq_base=`basename ${fq}`
	bzcat ${fq_base} > ${fq_base%.bz2}
	gzip -c ${fq_base%.bz2} > ${fq_base%.bz2}.gz
	rm $fq_base ${fq_base%.bz2}
done
#i=10
#while [ $i -le 10 ]; do
#cat 141023_HSGA.IMKB.IMKB${i}.fastq 141110_HSGA.IMKB.IMKB${i}.fastq > HSGA.IMKB.IMKB${i}.fastq
#gzip -c HSGA.IMKB.IMKB${i}.fastq > HSGA.IMKB.IMKB${i}.fastq.gz
#let i=i+1
#rm *.fastq.bz2 *.fastq
#done
