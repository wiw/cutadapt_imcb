#!/bin/bash
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #CTCGGCCGCGACC
i=12
adptr5=`echo $ADPTR_SHORT_5 | sed -e 's/^.\{'$i'\}//'` # Удаляю основания от адаптера с 5' конца
adptr3=`echo $ADPTR_SHORT_3 | sed -e 's/.\{'$i'\}$//'` # Удаляю основания от адаптера с 3' конца

cat $1 | perl -e 'while(<STDIN>) {$adptr5=$ARGV[0]; $adptr3=$ARGV[1]; $adptr5gatc=$adptr5.$ARGV[2]; $adptr3gatc=$adptr3.$ARGV[2]; $l=length($adptr5);  $id1=<STDIN>; $seq=<STDIN>; $id2=<STDIN>; $qual=<STDIN>; if($seq =~ /$adptr3gatc$/) {$seq=substr($seq, 1, -$l)."\n", length($seq); $qual=substr($qual,1, -$l) . "\n" }; if($seq =~ /^$adptr5gatc/){$seq=substr($seq,$l); $qual=substr($qual,$l)};  print $id1,$seq,$id2,$qual }' ${adptr5} ${adptr3} "GATC" > $2
