#!/bin/bash
dir=/home/anton/backup/cuted
fastq_files=$dir/*.fastq
cd $dir
for fq in $fastq_files; do
	{
	fq_base=`basename ${fq}`
	gzip -c ${fq_base} > ${fq_base}.gz
} &
done
wait
