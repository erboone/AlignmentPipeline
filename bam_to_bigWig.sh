#!/bin/bash

outdir=$2
input_path=$1

for fl in `ls $input_path/*_sorted.bam | sed 's|.*/||'`; do
name=${fl::-11}
echo "Converting: "$name" to bigWig file."
echo "Saving as: "$outdir"/"{$name}.bam"
bedtools genomeCoverageBed -d -bga -g ./ -ibam $input_path/${name}_sorted.bam  > $outdir/${name}.bam
done

