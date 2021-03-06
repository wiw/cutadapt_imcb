#!/bin/bash
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #CTCGGCCGCGACC
i=12
adptr5=`echo $ADPTR_SHORT_5 | sed -e 's/^.\{'$i'\}//'` # Удаляю основания от адаптера с 5' конца
adptr3=`echo $ADPTR_SHORT_3 | sed -e 's/.\{'$i'\}$//'` # Удаляю основания от адаптера с 3' конца

cat $1 | perl -e 'while(<STDIN>) { ($adptr5,$adptr3,$site)=@ARGV; $adptr5gatc=$adptr5.$site; $adptr3gatc=$site.$adptr3; $l=length($adptr5); $id1=$_; $seq=<STDIN>; $id2=<STDIN>; $qual=<STDIN>; if($seq =~ /$adptr3gatc$/) {$seq=substr($seq, 1, -$l-1)."\n"; $qual=substr($qual,1, -$l-1) . "\n" }; if($seq =~ /^$adptr5gatc/){$seq=substr($seq,$l); $qual=substr($qual,$l)}; print "$id1$seq$id2$qual" }' ${adptr5} ${adptr3} "GATC" > $2
