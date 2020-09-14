#!/usr/bin/bash

input=$1
base=$(basename $input .txt)

## filter low reads
awk '{if($4 > 1) print $0}' $input > ${base}.filtered

## fix start position (+1) compared to circexplorer2, find_circ, circRNA_finder
while read line 
        do

        chr=$(echo $line | awk '{print $1}')
        start=$(echo $line | awk '{sum = $2 -1; print sum}')
        stop=$(echo $line | awk '{print $3}')
        sign=$(echo $line | awk '{print $5}')
        count=$(echo $line | awk '{print $4}')
        
        echo -e "$chr\t$start\t$stop\t$sign\t$count" >> ${base}.bed

done < ${base}.filtered
