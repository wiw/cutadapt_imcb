#!/bin/bash
FASTQ_FILES=/home/anton/data/SRC/2047_3_1.3.1_Dam_Scrb__1.3.1_Dam_Scrb__ACAGTGA_L003_R1_001.fastq.gz
FASTX_REVCOM=fastx_reverse_complement
DIR=/home/anton/data/debug/temp
cd $DIR

# Make temporary files
out=`mktemp tmp.XXXXXXXXXX`
untrim_out=`mktemp tmp.XXXXXXXXXX`
out_wo_adapt_gatcs_len9=`mktemp tmp.XXXXXXXXXX`
untrim_out_gatcs_orig_len=`mktemp tmp.XXXXXXXXXX`
untrim_out_wo_gatcs_orig_len=`mktemp tmp.XXXXXXXXXX`

ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #CTCGGCCGCGACC
for fq in ${FASTQ_FILES}; do
	L=13
	fq_base=`basename ${fq}`
	#while [ ${L} -le 13 ]; do

	# Main trim reads
	cutadapt -g "${ADPTR_SHORT_5}" -a "${ADPTR_SHORT_3}" -n 3 -O ${L} --match-read-wildcards --untrimmed-output $DIR/${untrim_out} ${fq} -o $DIR/${out} >> $DIR/clip_${L}_${fq_base%.fastq.gz}.stats

	# Remove reads smaller then 9 bp
	cutadapt -b "GATC" -O 4 -m 9 --no-trim $DIR/${out} -o $DIR/${out_wo_adapt_gatcs_len9} >> $DIR/clip_${L}_${fq_base%.fastq.gz}.stats

	# Sort reads in untrimmed reads by presence GATC's
	cutadapt -b "GATC" -O 4 --no-trim --untrimmed-output $DIR/${untrim_out_wo_gatcs_orig_len} $DIR/${untrim_out} -o $DIR/${untrim_out_gatcs_orig_len} >> $DIR/clip_${L}_${fq_base%.fastq.gz}.stats

	#Sort reads by edge and with shift from 1 to 9 bases of adapters
	count=5
	p1=1
	p2=2
	cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $DIR/$p1/inner0-gatcs.fastq $DIR/${out_wo_adapt_gatcs_len9} -o $DIR/$p1/output0-gatcs.fastq > $DIR/$p1/clip_${out_wo_adapt_gatcs_len9}_gatcs.stats
	
	cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $DIR/$p2/inner0-gatcs.fastq $DIR/${untrim_out_gatcs_orig_len} -o $DIR/$p2/output0-gatcs.fastq > $DIR/$p2/clip_${untrim_out_gatcs_orig_len}_gatcs.stats
		while [ $count -le 13 ]; do
			i=$((13-$((${count}-4)))) # start from 12
			adptr5=`echo $ADPTR_SHORT_5 | sed -e 's/^.\{'$i'\}//'`
			adptr3=`echo $ADPTR_SHORT_3 | sed -e 's/.\{'$i'\}$//'`
						
			# search gatcs from reads w/o adapters, with gatcs, length >= 9
			cutadapt -g "^${adptr5}GATC" -a "GATC${adptr3}$" -O $count -e 0.01 --no-trim --untrimmed-output $DIR/$p1/inner$((${count}-4))-gatcs.fastq $DIR/$p1/inner$((${count}-5))-gatcs.fastq -o $DIR/$p1/output$((${count}-4))-gatcs.fastq >> clip_${out_wo_adapt_gatcs_len9}_gatcs.stats

			# search gatcs from reads with gatcs, original length
			cutadapt -g "^${adptr5}GATC" -a "GATC${adptr3}$" -O $count -e 0.01 --no-trim --untrimmed-output $DIR/$p2/inner$((${count}-4))-gatcs.fastq $DIR/$p2/inner$((${count}-5))-gatcs.fastq -o $DIR/$p2/output$((${count}-4))-gatcs.fastq >> clip_${out_wo_adapt_gatcs_len9}_gatcs.stats

count=$(($count+1)) # 6
	done
	#grep "Trimmed reads" clip_${L}_${fq_base%.fastq.gz}.stats | sed -e 's/^..................................\(\)/\1/' | sed 's/..$//' >> clip_common_${fq_base%.fastq.gz}.stats
		#rm clip_${L}_${fq_base%.fastq.gz}.stats
		#L=$(($L+1))
	#done
done
