#!/usr/bin/bash 

input=$1
base=$(basename $input .gtf)

## strip headers from file and remove BSJ with less than 2 reads
grep -v "#" $input | grep -v "bsj 1.000" > ${base}.filtered
awk '{print \$14}' ${base}.filtered | cut -d'.' -f1 > counts
cut -f 1,4,5,7 ${base}.filtered > ${base}.txt
paste ${base}.txt counts > ${base}.bed.tmp

## fix start position (+1) compared to circexplorer2, find_circ, circRNA_finder
while read line 
        do

        chr=$(echo $line | awk '{print $1}')
        start=$(echo $line | awk '{sum = $2 -1; print sum}')
        stop=$(echo $line | awk '{print $3}')
        sign=$(echo $line | awk '{print $4}')

        echo -e "$chr\t$start\t$stop\t$sign" >> ${base}.bed

done < ${base}.bed.tmp
