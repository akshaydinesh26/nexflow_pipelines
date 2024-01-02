#!/usr/bin/bash

read1=$1;
read2=$2;
outname=$3;
trimval=$4;
len=$5;
cores=$6;

fastp -i ${reads1} -I ${reads2} -o "filt_${outname}_1.fastq.gz" -O "filt_${outname}_2.fastq.gz" -w $cores -q $trimval -l $len -R "${outname}" -j "${outname}.json" -h "${outname}.html"