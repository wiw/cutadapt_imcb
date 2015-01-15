#!/bin/bash
FASTQ_FILES=/home/anton/backup/input/1540*.fastq.gz # расположение исходных файлов
FASTX_REVCOM=fastx_reverse_complement 
DIR=/home/anton/backup/output # папка куда будет всё складываться
#OUT=/home/anton/backup/output # папка куда будут перемещаться данные
DSCR=/home/anton/data/damid_description.csv # расположение файла описаний

# Определяем адаптеры
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #CTCGGCCGCGACC
ILLUMINA_5="GCTCTTCCGATCT"
ILLUMINA_3=`echo -e ">\n${ILLUMINA_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #AGATCGGAAGAGC

# Для каждого файла из FASTQ_FILES выполняется следующая конструкция
for fq in ${FASTQ_FILES}; do
{
	fq_base=`basename ${fq}` # Оставляем только имя файла
		fq_human=`grep -w $fq_base $DSCR | sed "s/[a-zA-Z0-9_-.]*$//;s/	//"` # человекопонятное имя файла, берем из файла описаний
		mkdir $DIR/${fq_base%.fastq.gz} # Создаем папку с именем файла, но без расширения
		basef=$DIR/${fq_base%.fastq.gz}
	stats=$DIR/${fq_base%.fastq.gz}/stats
		len9=$DIR/${fq_base%.fastq.gz}/len9 # Имя подпапки куда будет складываться весь результат для обрезанных ридов
		olen=$DIR/${fq_base%.fastq.gz}/orig_len # Имя подпапки куда будет складываться весь результат для ридов оригинальной длины
		mkdir $len9
		mkdir $olen
		mkdir $stats
		count=5 # Соответствует длине адаптера с GATC фрагментом +1 основание (итого 5)

# Main trim reads. Обработка cutadapt со следующими параметрами: 5' & 3' адаптеры встречаются от 3-х и более раз, перекрытие от 9-ти и более оснований. Те риды где найдены адаптеры уходят в файл out.fastq, где не нашлись адаптеры - в файл untrim_out.fastq
		cutadapt -g "${ADPTR_SHORT_5}" -g "${ILLUMINA_5}" -a "${ADPTR_SHORT_3}" -a "${ILLUMINA_3}" -n 3 -O 9 --match-read-wildcards --untrimmed-output $basef/untrim_out.fastq ${fq} -o $basef/out.fastq > $stats/clip_${fq_base%.fastq.gz}.stats

# Remove reads smaller then 9 bp. Обрабатываю файл с отрезанными адаптерами - ищу GATC фрагменты в любом положении по ридам с минимальной длиной 9 оснований, ничего не отрезаю. Что нашлось уходит в файл out_wo_adapt_gatcs_len9.fastq, риды с меньшей длиной и/или без GATC уходят в мусорку
		cutadapt -g "GATC" -a "GATC" -O 4 -m 9 --no-trim --untrimmed-output $basef/out_wo_adapt_wo_gatcs_small_len.fastq $basef/out.fastq -o $basef/out_wo_adapt_gatcs_len9.fastq > $stats/clip_len9_${fq_base%.fastq.gz}.stats

# Sort reads in untrimmed reads by presence GATC's. Обрабатываю файл с ридами без адаптеров - также ищу GATC фрагменты. Ничего не отрезаю. Риды с фрагментами отправляются в файл untrim_out_gatcs_orig_len.fastq, без фрагментов идут в файл untrim_out_wo_gatcs_orig_len.fastq.
		cutadapt -g "GATC" -a "GATC" -O 4 --no-trim --untrimmed-output $basef/untrim_out_wo_gatcs_orig_len.fastq $basef/untrim_out.fastq -o $basef/untrim_out_gatcs_orig_len.fastq > $stats/clip_orig_len_${fq_base%.fastq.gz}.stats

#Sort reads by edge and with shift from 1 to 9 bases of adapters. Этот блок кода вычленяет риды с GATC фрагментами на краях и на небольшом отступе от края - вплоть до 9-ти нуклеотидов из адаптера, как с прямой так и с обратной стороны.

# Первоначальная обработка и поиск только крайних фрагментов из файла out_wo_adapt_gatcs_len9.fastq, перекрытие 4 основания, 1% ошибок
		cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $len9/inner0-gatcs.fastq $basef/out_wo_adapt_gatcs_len9.fastq -o $len9/output0-gatcs.fastq > $stats/clip_len9_gatcs4.stats

# Первоначальная обработка и поиск только крайних фрагментов из файла untrim_out_gatcs_orig_len.fastq, перекрытие 4 основания, 1% ошибок
		cutadapt -g "^GATC" -a "GATC$" -O 4 -e 0.01 --no-trim --untrimmed-output $olen/inner0-gatcs.fastq $basef/untrim_out_gatcs_orig_len.fastq -o $olen/output0-gatcs.fastq > $stats/clip_orig_len_gatcs4.stats

#############################
### Переменные для отчета ###
#############################
		s0_reads=`grep "Processed reads" $stats/clip_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//'`

		s1_trim=`grep "Trimmed reads" $stats/clip_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`
		s1_untrim=$((${s0_reads} - ${s1_trim}))

		s2_trim_gatcs=$((`wc -l $basef/out_wo_adapt_gatcs_len9.fastq | sed 's/[a-zA-Z0-9_./-]*$//'`/4))
		s2_untrim_trash_gatcs=$((`wc -l $basef/out_wo_adapt_wo_gatcs_small_len.fastq | sed 's/[a-zA-Z0-9_./-]*$//'`/4)) #!!!!!!!
#s2_rem_too_short=`grep "Too short reads" $stats/clip_len9_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%)0-9.a-z ]*$//;s/[( ]*$//'`
#s2_rem_match_reads=`grep "Matched reads" $stats/clip_len9_${fq_base%.fastq.gz}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`
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
#s2_rem_too_short_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_rem_too_short; (b/a)*100" | sed 's/[0].$//'`%
#s2_rem_match_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_rem_match_reads; (b/a)*100" | sed 's/[0].$//'`%
		s2_trash_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_trash_reads; (b/a)*100" | sed 's/[0].$//'`%
		s2_untrim_gatc_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_untrim_gatc; (b/a)*100" | sed 's/[0].$//'`%
		s2_untrim_pct=`bc <<< "scale=4; a=$s0_reads; b=$s2_untrim; (b/a)*100" | sed 's/[0].$//'`%

		s3_input_trim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_input_trim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_match_trim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_match_trim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_input_untrim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_input_untrim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_match_untrim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_match_untrim_reads; (b/a)*100" | sed 's/[0].$//'`%

		echo "<tr><td>GATC</td><td><script>document.write(number_format(${s3_input_trim_reads}, 0, '.', ' '))</script> (${s3_input_trim_reads_pct})</td><td><script>document.write(number_format(${s3_match_trim_reads}, 0, '.', ' '))</script> (${s3_match_trim_reads_pct})</td></tr>" >> $stats/len9.txt
		echo "<tr><td>GATC</td><td><script>document.write(number_format(${s3_input_untrim_reads}, 0, '.', ' '))</script> (${s3_input_untrim_reads_pct})</td><td><script>document.write(number_format(${s3_match_untrim_reads}, 0, '.', ' '))</script> (${s3_match_untrim_reads_pct})</td></tr>" >> $stats/orig_len.txt

# Запуск цикла в котором идет поиск фрагментов GATC отстоящих от края рида на 1, 2, 3 и так далее до 9-ти оснований
		while [ $count -le 13 ]; do # 4 + 9 = 13 длина короткого адаптера 
			i=$((13-$((${count}-4)))) # start from 12 Считаю число оснований, которые нужно удалить от оригинального адаптера
				adptr5=`echo $ADPTR_SHORT_5 | sed -e 's/^.\{'$i'\}//'` # Удаляю основания от адаптера с 5' конца
	adptr3=`echo $ADPTR_SHORT_3 | sed -e 's/.\{'$i'\}$//'` # Удаляю основания от адаптера с 3' конца

# search gatcs from reads w/o adapters, with gatcs, length >= 9. Поиск отстоящих GATC фрагментов на ридах без адаптеров, длиной не менее 9 оснований
		cutadapt -g "^${adptr5}GATC" -a "GATC${adptr3}$" -O $count -e 0.01 --no-trim --untrimmed-output $len9/inner$((${count}-4))-gatcs.fastq $len9/inner$((${count}-5))-gatcs.fastq -o $len9/output$((${count}-4))-gatcs.fastq > $stats/clip_len9_gatcs${count}.stats

# Удаление огрызков адаптера для найденых ридов и вывод в отдельных файл 
		cat $len9/output$((${count}-4))-gatcs.fastq | sed -e "s/^"$adptr5"GATC/GATC/" | sed -e "s/GATC"${adptr3}"$/GATC/" > $len9/sed_output${count}-gatcs.fastq

# search gatcs from reads with gatcs, original length. Поиск отстоящих GATC фрагментов на ридах без адаптеров, оригинальной длины
		cutadapt -g "^${adptr5}GATC" -a "GATC${adptr3}$" -O $count -e 0.01 --no-trim --untrimmed-output $olen/inner$((${count}-4))-gatcs.fastq $olen/inner$((${count}-5))-gatcs.fastq -o $olen/output$((${count}-4))-gatcs.fastq > $stats/clip_orig_len_gatcs${count}.stats

# Удаление огрызков адаптера для найденых ридов и вывод в отдельных файл 
		cat $olen/output$((${count}-4))-gatcs.fastq | sed -e "s/^"$adptr5"GATC/GATC/" | sed -e "s/GATC"${adptr3}"$/GATC/" > $olen/sed_output${count}-gatcs.fastq

		s3_input_trim_reads=`grep "Processed reads" $stats/clip_len9_gatcs${count}.stats | sed 's/^[a-zA-Z ^t:]*//'`
		s3_match_trim_reads=`grep "Matched reads" $stats/clip_len9_gatcs${count}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`
		s3_input_untrim_reads=`grep "Processed reads" $stats/clip_orig_len_gatcs${count}.stats | sed 's/^[a-zA-Z ^t:]*//'`
		s3_match_untrim_reads=`grep "Matched reads" $stats/clip_orig_len_gatcs${count}.stats | sed 's/^[a-zA-Z ^t:]*//;s/[%()0-9.]*$//;s/[ ^]*$//'`

		s3_input_trim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_input_trim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_match_trim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_match_trim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_input_untrim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_input_untrim_reads; (b/a)*100" | sed 's/[0].$//'`%
		s3_match_untrim_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s3_match_untrim_reads; (b/a)*100" | sed 's/[0].$//'`%

		echo "<tr><td>${adptr5}GATC</td><td><script>document.write(number_format(${s3_input_trim_reads}, 0, '.', ' '))</script> (${s3_input_trim_reads_pct})</td><td><script>document.write(number_format(${s3_match_trim_reads}, 0, '.', ' '))</script> (${s3_match_trim_reads_pct})</td></tr>" >> $stats/len9.txt
		echo "<tr><td>${adptr5}GATC</td><td><script>document.write(number_format(${s3_input_untrim_reads}, 0, '.', ' '))</script> (${s3_input_untrim_reads_pct})</td><td><script>document.write(number_format(${s3_match_untrim_reads}, 0, '.', ' '))</script> (${s3_match_untrim_reads_pct})</td></tr>" >> $stats/orig_len.txt

		count=$(($count+1)) # 6 Увеличиваю счетчик на единицу
		done
# Объединение всех найденых ридов в один файл
		cat $len9/output0-gatcs.fastq $len9/sed_output*-gatcs.fastq $olen/output0-gatcs.fastq $olen/sed_output*-gatcs.fastq > $basef/interim_gatcs_${fq_base}.fastq
#cat $olen/output0-gatcs.fastq $olen/sed_output*-gatcs.fastq > $olen/untrim_${fq_base}.fastq

# Удаление внутренних ридов в файле с краевыми GATC's
		pre=`head -n 1 $basef/interim_gatcs_${fq_base}.fastq | cut -c 1-2`
		grep -E -v '.+(GATC)+.+' $basef/interim_gatcs_${fq_base}.fastq | sed '/^[-].*$/d' | sed ":loop; N; /\n+/ ! { $ ! b loop }; /\n${pre}[^\n]\+\n+/ d" > $basef/summary_gatcs_${fq_base}.fastq

		s4_interim_gatcs=$((`wc -l $basef/interim_gatcs_${fq_base}.fastq | sed 's/[a-zA-Z0-9_./-]*$//'`/4))
		s4_interim_gatcs_pct=`bc <<< "scale=4; a=$s0_reads; b=$s4_interim_gatcs; (b/a)*100" | sed 's/[0].$//'`%

		s4_interim_trash_reads=$((${s3_input_trim_reads}+${s3_input_untrim_reads}+${s2_trash_reads}-${s3_match_trim_reads}-${s3_match_untrim_reads}))
		s4_interim_trash_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s4_interim_trash_reads; (b/a)*100" | sed 's/[0].$//'`%

		s3_html_trim=`cat $stats/len9.txt`
		s3_html_untrim=`cat $stats/orig_len.txt`

		s5_summary_gatcs=$((`wc -l $basef/summary_gatcs_${fq_base}.fastq | sed 's/[a-zA-Z0-9_./-]*$//'`/4))
		s5_summary_gatcs_pct=`bc <<< "scale=4; a=$s0_reads; b=$s5_summary_gatcs; (b/a)*100" | sed 's/[0].$//'`%

		s5_trash_reads=$((${s4_interim_trash_reads}+${s5_summary_gatcs}-${s4_interim_gatcs}))
		s5_trash_reads_pct=`bc <<< "scale=4; a=$s0_reads; b=$s5_trash_reads; (b/a)*100" | sed 's/[0].$//'`%

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
		<p align=\"center\">original length<\br>without GATC</p>
		</div>
		<div class=\"span4\" style=\"background-color: #448f30\">
		<h4 align=\"center\"><script>document.write(number_format(${s4_interim_gatcs}, 0, '.', ' '))</script> edge reads (${s4_interim_gatcs_pct})</h4>
		<p align=\"center\">original length with edge GATC(s) & truncated<\br>with GATC(s) at the edge(s)</p>
		</div>
		<div class=\"span4\" style=\"color: #f10026\">
		<h4 align=\"center\"><script>document.write(number_format(${s4_interim_trash_reads}, 0, '.', ' '))</script> trash reads (${s4_interim_trash_reads_pct})</h4>
		<p align=\"center\">too short after removal of adapters<\br>no GATC(s) after removal of adapters<\br>contain inner GATS(s)</p>
		</div>
		</div>
		<p>&nbsp;</p>
		<div class=\"row\">
		<div class=\"alert\">Removing  reads with inner GATC(s)</div>
		<div class=\"span4\">
		<h4 align=\"center\"><script>document.write(number_format(${s2_untrim}, 0, '.', ' '))</script> inner reads (${s2_untrim_pct})</h4>
		<p align=\"center\">original length<\br>without GATC</p>
		</div>
		<div class=\"span4\" style=\"background-color: #448f30\">
		<h4 align=\"center\"><script>document.write(number_format(${s5_summary_gatcs}, 0, '.', ' '))</script> edge reads (${s5_summary_gatcs_pct})</h4>
		<p align=\"center\">original length with edge GATC(s) & truncated<\br>with GATC(s) at the edge(s)</p>
		</div>
		<div class=\"span4\" style=\"color: #f10026\">
		<h4 align=\"center\"><script>document.write(number_format(${s5_trash_reads}, 0, '.', ' '))</script> trash reads (${s5_trash_reads_pct})</h4>
		<p align=\"center\">too short after removal of adapters<\br>no GATC(s) after removal of adapters<\br>contain inner GATS(s)</p>
		</div>
		</div>

		<hr>
		<footer><b><p>&copy; Laboratory of cell division IMCB SB RAS, Novosibirsk, 2014-2015</p></b></footer>
		</div>
		</body>
		</html>" > $stats/${fq_human}_report.html


		rm $len9/inner*.fastq $len9/output*.fastq $len9/sed_*.fastq $olen/inner*.fastq $olen/output*.fastq $olen/sed_*.fastq $basef/out*.fastq $basef/untrim_out.fastq $basef/untrim_out_gatcs_orig_len.fastq $basef/interim_gatcs_${fq_base}.fastq
#mv $basef $OUT
} &
done
wait

