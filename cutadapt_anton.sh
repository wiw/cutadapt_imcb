#!/bin/bash
FASTQ_FILES=$1/*.fastq.gz # Source data in fastq.gz format
FASTX_REVCOM=fastx_reverse_complement # path to fastx_reverse_complement
DIR=$2 # Output dir
#OUT=/home/anton/backup/output # If we have fast SSD with big volume space then output dir may be folder to SSD, and $OUT folder may be locate on big HDD. Data at first are recorded on SSD then move to HDD. This option disabled by default.
DSCR=/home/anton/data/DAM/RUN/damid_description.csv # path to description file which establishes a correspondence between the file name and its human-readable name.

# Define adapters
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #CTCGGCCGCGACC
ILLUMINA_5="GCTCTTCCGATCT"
ILLUMINA_3=`echo -e ">\n${ILLUMINA_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #AGATCGGAAGAGC

# Start statistic record
echo "Data.set;fastq.file;Total.number.of.reads.obtained;Without.GATCs.original.length;With.edge.GATCs;Cutadapt.Trash;Sum;Sum-Total.number.of.reads.obtained;Without.GATCs.original.length.Percentage;With.edge.GATCs.Percentage;Cutadapt.Trash.Percentage" > $DIR/statistics.csv

# For each file in FASTQ_FILES run this block 
for fq in ${FASTQ_FILES}; do
{
	fq_base=`basename ${fq}` # save only file name
		fq_human=`grep -w $fq_base $DSCR | sed "s/[a-zA-Z0-9_-.]*$//;s/	//"` # human-readable name in variable, get by parse $DSCR
		mkdir $DIR/${fq_base%.fastq.gz} # make dir named as filename
		basef=$DIR/${fq_base%.fastq.gz}
		stats=$DIR/${fq_base%.fastq.gz}/stats
		len9=$DIR/${fq_base%.fastq.gz}/len9 # folder for cuted reads 
		olen=$DIR/${fq_base%.fastq.gz}/orig_len # folder for uncuted reads
# make folder defined on last step
		mkdir $len9
		mkdir $olen
		mkdir $stats
		count=5 # some variables corresponds to the length of the adapter with a fragment of one GATC base (5 total)

# Main trim reads. Processing cutadapt with the following parameters: 5 '& 3' adapters encountered from 3 or more times, the overlap of the adapter 9 or more bases. Those reads where found adapters write to file out.fastq, where there was no adapters write to file untrim_out.fastq 
		cutadapt -g "${ADPTR_SHORT_5}" -g "${ILLUMINA_5}" -a "${ADPTR_SHORT_3}" -a "${ILLUMINA_3}" -e 0.01 -n 3 --overlap 12 --match-read-wildcards --untrimmed-output $basef/untrim_out.fastq --too-short $basef/s1-too-short.fastq ${fq} -o $basef/out.fastq > $stats/clip_${fq_base%.fastq.gz}.stats

# Remove reads smaller then 9 bp. Process files with truncated adapters - looking GATC fragments in any position by reads with a minimum length of 9 bases, do not cut off. That there was a goes into file out_wo_adapt_gatcs_len9.fastq, reads with smaller length and / or without GATC goes in the trash 
		cutadapt -g "GATC" -a "GATC" -O 4 -m 9 --no-trim --untrimmed-output $basef/out_wo_adapt_wo_gatcs_small_len.fastq --too-short $basef/s2-too-short.fastq $basef/out.fastq -o $basef/out_wo_adapt_gatcs_len9.fastq > $stats/clip_len9_${fq_base%.fastq.gz}.stats

# Sort reads in untrimmed reads by presence GATC's. Process files with reads without adapters - am also looking for GATC fragments. Nothing is cut off. Reads with fragments are sent to a file untrim_out_gatcs_orig_len.fastq, without going into the file fragments untrim_out_wo_gatcs_orig_len.fastq. 
		cutadapt -g "GATC" -a "GATC" -O 4 --no-trim --untrimmed-output $basef/untrim_out_wo_gatcs_orig_len.fastq $basef/untrim_out.fastq -o $basef/untrim_out_gatcs_orig_len.fastq > $stats/clip_orig_len_${fq_base%.fastq.gz}.stats

# Sort reads by edge and with shift from 1 to 9 bases of adapters. This block of code looking for reads with GATC fragments at the edges and slightly indented from the edge - up to 9 nucleotides from the adapter as a direct and the reverse side. 

# Initial processing and search only the edge of the file fragments out_wo_adapt_gatcs_len9.fastq, overlapping 4 base, 1% error
		cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $len9/inner0-gatcs.fastq $basef/out_wo_adapt_gatcs_len9.fastq -o $len9/output0-gatcs.fastq > $stats/clip_len9_gatcs4.stats

# Initial processing and search only the edge of the file fragments untrim_out_gatcs_orig_len.fastq, overlapping 4 base 1% error
		cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $olen/inner0-gatcs.fastq $basef/untrim_out_gatcs_orig_len.fastq -o $olen/output0-gatcs.fastq > $stats/clip_orig_len_gatcs4.stats

#############################
###  Variable for report  ###
#############################
		s0_reads=`grep "Processed reads" $stats/clip_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//'`

		s1_trim=`grep "Trimmed reads" $stats/clip_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`
		s1_untrim=$((${s0_reads} - ${s1_trim}))

		s2_trim_gatcs=`grep "^\+$" $basef/out_wo_adapt_gatcs_len9.fastq | wc -l`
		s2_untrim_trash_gatcs=`grep "^\+$" $basef/out_wo_adapt_wo_gatcs_small_len.fastq | wc -l`
		s2_trash_reads=$(($s1_trim-$s2_trim_gatcs))
		s2_untrim_gatc=`grep "Matched reads" $stats/clip_orig_len_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`
		s2_untrim=$(($s1_untrim-$s2_untrim_gatc))

		s3_input_trim_reads=`grep "Processed reads" $stats/clip_len9_gatcs4.stats | sed 's/^[a-zA-Z ^t:]*//'`
		s3_match_trim_reads=`grep "Matched reads" $stats/clip_len9_gatcs4.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`
		s3_input_untrim_reads=`grep "Processed reads" $stats/clip_orig_len_gatcs4.stats | sed 's/^[a-zA-Z ^t:]*//'`
		s3_match_untrim_reads=`grep "Matched reads" $stats/clip_orig_len_gatcs4.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`

# Calculate percentages
		s0_reads_pct=`bc <<< "a=$s0_reads; b=$s0_reads; (b/a)*100"`%
		s1_trim_pct=`bc <<< "scale=4; a=$s0_reads; b=$s1_trim; (b/a)*100" | sed 's/[0].$//'`%
		s1_untrim_pct=`bc <<< "scale=4; a=$s0_reads; b=$s1_untrim; (b/a)*100" | sed 's/[0].$//'`%
		s2_trim_gatcs_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_trim_gatcs; (b/a)*100" | sed 's/[0].$//'`%
		s2_untrim_trash_gatcs_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_untrim_trash_gatcs; (b/a)*100" | sed 's/[0].$//'`%
		s2_trash_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_trash_reads; (b/a)*100" | sed 's/[0].$//'`%
		s2_untrim_gatc_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_untrim_gatc; (b/a)*100" | sed 's/[0].$//'`%
		s2_untrim_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_untrim; (b/a)*100" | sed 's/[0].$//'`%

		s3_input_trim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_input_trim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_match_trim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_match_trim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_input_untrim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_input_untrim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_match_untrim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_match_untrim_reads; (b/a)*100" | sed 's/[0].$//'`%
#############################

# record statistic to temporary file
		echo "<tr><td>GATC</td><td><script>document.write(number_format(${s3_input_trim_reads}, 0, '.', ' '))</script> (${s3_input_trim_reads_pct})</td><td><script>document.write(number_format(${s3_match_trim_reads}, 0, '.', ' '))</script> (${s3_match_trim_reads_pct})</td></tr>" >> $stats/len9.txt
		echo "<tr><td>GATC</td><td><script>document.write(number_format(${s3_input_untrim_reads}, 0, '.', ' '))</script> (${s3_input_untrim_reads_pct})</td><td><script>document.write(number_format(${s3_match_untrim_reads}, 0, '.', ' '))</script> (${s3_match_untrim_reads_pct})</td></tr>" >> $stats/orig_len.txt


# Combine all founded reads to one file
		cat $len9/output0-gatcs.fastq $olen/output0-gatcs.fastq > $basef/interim_gatcs_${fq_base}.fastq

# Remove reads with inner GATC's
		pre=`head -n 1 $basef/interim_gatcs_${fq_base}.fastq | cut -c 1-2`
		sed -r "s/.+(GATC)+.+/empty sequence/" $basef/interim_gatcs_${fq_base}.fastq | sed "/^$pre/ {N; /empty sequence/ { N; /\n+$/ { N; d } } }" > $basef/summary_gatcs_${fq_base}.fastq

#############################
###  Variable for report  ###
#############################
		s4_interim_gatcs=`grep "^\+$" $basef/interim_gatcs_${fq_base}.fastq | wc -l`
		s4_interim_gatcs_pct=`bc <<< "scale=4; a=$s0_reads; b=$s4_interim_gatcs; (b/a)*100" | sed 's/[0].$//'`%

		s4_interim_trash_reads=$((${s0_reads}-${s4_interim_gatcs}-${s2_untrim}))
		s4_interim_trash_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s4_interim_trash_reads; (b/a)*100" | sed 's/[0].$//'`%

		s3_html_trim=`cat $stats/len9.txt`
		s3_html_untrim=`cat $stats/orig_len.txt`

		s5_summary_gatcs=`grep "^\+$" $basef/summary_gatcs_${fq_base}.fastq | wc -l`
		s5_summary_gatcs_pct=`bc <<< "scale=4; a=$s0_reads; b=$s5_summary_gatcs; (b/a)*100" | sed 's/[0].$//'`%

		s5_trash_reads=$((${s0_reads}-${s5_summary_gatcs}-${s2_untrim}))
		s5_trash_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s5_trash_reads; (b/a)*100" | sed 's/[0].$//'`%
#############################

# Make text statistic in csv file
echo "$fq_human;$fq_base;$s0_reads;$s2_untrim;$s5_summary_gatcs;$s5_trash_reads;$(($s2_untrim+$s5_summary_gatcs+$s5_trash_reads));$(($s0_reads-$(($s2_untrim+$s5_summary_gatcs+$s5_trash_reads))));${s2_untrim_pct%\%};${s5_summary_gatcs_pct%\%};${s5_trash_reads_pct%\%}" >> $DIR/statistics.csv

# Make visual statistic in html file
		echo "
		<!DOCTYPE html>
		<html lang=\"en\">
		<head>
		<meta charset=\"utf-8\">
		<title>Report on the work cutadapt program</title>
		<style type=\"text/css\">body{padding-top: 60px; padding-bottom: 40px;}</style>
		<link href=\"http://cellbiol.ru/files/bootstrap/css/bootstrap.css\" rel=\"stylesheet\">
		<link href=\"http://cellbiol.ru/files/bootstrap/css/bootstrap-responsive.css\" rel=\"stylesheet\">
		<script>
		function number_format( number, decimals, dec_point, thousands_sep ) {	// Format a number with grouped thousands
			// 
			// +   original by: Jonas Raoni Soares Silva (http://www.jsfromhell.com)
			// +   improved by: Kevin van Zonneveld (http://kevin.vanzonneveld.net)
			// +	 bugfix by: Michael White (http://crestidg.com)

			var i, j, kw, kd, km;

			// input sanitation & defaults
			if( isNaN(decimals = Math.abs(decimals)) ){
				decimals = 2;
			}
			if( dec_point == undefined ){
				dec_point = \",\";
			}
			if( thousands_sep == undefined ){
				thousands_sep = \".\";
			}

			i = parseInt(number = (+number || 0).toFixed(decimals)) + \"\";

			if( (j = i.length) > 3 ){
				j = j % 3;
			} else{
				j = 0;
			}

			km = (j ? i.substr(0, j) + thousands_sep : \"\");
			kw = i.substr(j).replace(/(\d{3})(?=\d)/g, \"\$1\" + thousands_sep);
			//kd = (decimals ? dec_point + Math.abs(number - i).toFixed(decimals).slice(2) : \"\");
			kd = (decimals ? dec_point + Math.abs(number - i).toFixed(decimals).replace(/-/, 0).slice(2) : \"\");


			return km + kw + kd;
		}</script>
	</head>
		<body>
		<div class=\"container\">
		<div class=\"hero-unit\">
		<h1 align=\"center\">Cutadapt report </br><small>from ${fq_base}</small></h1>
		</div>
		<div class=\"row\"> 
		<div class=\"span12\"> 
		<h2 align=\"center\">input: <script>document.write(number_format(${s0_reads}, 0, '.', ' '))</script> reads (${s0_reads_pct})</h2> 
		</div>
		</div>
		<div class=\"row\">
		<div class=\"alert\">Removing DamID and Illumina adapters</div>
		<div class=\"span6\"> <h3 align=\"center\"><script>document.write(number_format(${s1_untrim}, 0, '.', ' '))</script> reads (${s1_untrim_pct})</h3> 
		<p align=\"center\">No adapters &ge;9 bp found</p>
		</div>
		<div class=\"span6\"> <h3 align=\"center\"><script>document.write(number_format(${s1_trim}, 0, '.', ' '))</script> reads (${s1_trim_pct})</h3>
		<p align=\"center\">Adapters &ge;9 bp found and removed</p>
		</div>
		</div>
		<div class=\"row\">
		<div class=\"alert\">Looking for GATC motives</div>
		<div class=\"span3\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_untrim}, 0, '.', ' '))</script> reads (${s2_untrim_pct})</h4>
		<p align=\"center\">without GATC(s)</p>
		</div>
		<div class=\"span3\" style=\"background-color: #83a136\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_untrim_gatc}, 0, '.', ' '))</script> reads (${s2_untrim_gatc_pct})</h4>
		<p align=\"center\">with GATC(s)</p>
		</div>
		<div class=\"span3\" style=\"background-color: #abec00\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_trim_gatcs}, 0, '.', ' '))</script> reads (${s2_trim_gatcs_pct})</h4>
		<p align=\"center\">with GATC(s)</p>
		</div>
		<div class=\"span3\" style=\"color: #f10026\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_trash_reads}, 0, '.', ' '))</script> trash reads (${s2_trash_reads_pct})</h4>
		<p align=\"center\">too short (&lt;9 bp) reads or/and do not contain GATC</p>
		</div>
		</div>
		<div class=\"row\">
		<div class=\"alert\">Determing location of GATC motives</div>
		<div class=\"span6\" style=\"background-color: #83a136\">
		<table class=\"table table-striped\" align=\"center\">
		<thead>
		<tr><th>Location</th><th>Processed reads</th><th>Matched reads</th></tr></thead>${s3_html_untrim}</table>
		</div>
		<div class=\"span6\" style=\"background-color: #abec00\">
		<table class=\"table table-striped\" align=\"center\">
		<thead><tr><th>Location</th><th>Processed reads</th><th>Matched reads</th></tr></thead>${s3_html_trim}</table>
		</div>
		</div>
		<p>&nbsp;</p>
		<div class=\"row\">
		<div class=\"alert\">Interim report</div>
		<div class=\"span4\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_untrim}, 0, '.', ' '))</script> inner reads (${s2_untrim_pct})</h4>
		<p align=\"center\"><ul><li>original length</li><li>without GATC</li></ul></p>
		</div>
		<div class=\"span4\" style=\"background-color: #448f30\">
		<h4 align=\"center\"><script>document.write(number_format(${s4_interim_gatcs}, 0, '.', ' '))</script> edge reads (${s4_interim_gatcs_pct})</h4>
		<p align=\"center\"><ul><li>original length with edge GATC(s) & truncated</li><li>with GATC(s) at the edge(s)</li></ul></p>
		</div>
		<div class=\"span4\" style=\"color: #f10026\">
		<h4 align=\"center\"><script>document.write(number_format(${s4_interim_trash_reads}, 0, '.', ' '))</script> trash reads (${s4_interim_trash_reads_pct})</h4>
		<p align=\"center\"><ul><li>too short after removal of adapters</li><li>no GATC(s) after removal of adapters</li><li>contain inner GATC(s)</li></ul></p>
		</div>
		</div>
		<p>&nbsp;</p>
		<div class=\"row\">
		<div class=\"alert\">Removing  reads with inner GATC(s)</div>
		<div class=\"span4\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_untrim}, 0, '.', ' '))</script> inner reads (${s2_untrim_pct})</h4>
		<p align=\"center\"><ul><li>original length</li><li>without GATC</li></ul></p>
		</div>
		<div class=\"span4\" style=\"background-color: #448f30\">
		<h4 align=\"center\"><script>document.write(number_format(${s5_summary_gatcs}, 0, '.', ' '))</script> edge reads (${s5_summary_gatcs_pct})</h4>
		<p align=\"center\"><ul><li>original length with edge GATC(s) & truncated</li><li>with GATC(s) at the edge(s)</li></ul></p>
		</div>
		<div class=\"span4\" style=\"color: #f10026\">
		<h4 align=\"center\"><script>document.write(number_format(${s5_trash_reads}, 0, '.', ' '))</script> trash reads (${s5_trash_reads_pct})</h4>
		<p align=\"center\"><ul><li>too short after removal of adapters</li><li>no GATC(s) after removal of adapters</li><li>contain inner GATC(s)</li></ul></p>
		</div>
		</div>

		<hr>
		<footer><b><p>&copy; Laboratory of cell division IMCB SB RAS, Novosibirsk, 2014-2015</p></b></footer>
		</div>
		</body>
		</html>" > $basef/${fq_human}_report.html

# remove intermediate files
#		rm -R $len9/* $olen/* $stats/* $basef/out*.fastq $basef/untrim_out.fastq $basef/untrim_out_gatcs_orig_len.fastq $basef/interim_gatcs_${fq_base}.fastq
#mv $basef $OUT # If folder $OUT is defined then to move output data from $DIR to $OUT
} & # all files are processed in parallel processes
done
wait
