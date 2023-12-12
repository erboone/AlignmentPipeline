#!/bin/bash

outdir=$2
input_path=$1

for fl in `ls $input_path/MT*/*.sam | sed 's|.*/||'`; do
name=${fl::-16}
echo "Converting: "$name" to BAM file."
echo "Saving as: "$outdir"/"{$name}.bam
samtools view -bS $input_path/$name/${name}_Aligned.out.sam  > $outdir/${name}.bam
samtools sort $outdir/${name}.bam -o $outdir/${name}_sorted.bam
done

