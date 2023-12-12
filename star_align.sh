#!/bin/bash

outdir=$2
input_path=$1

for fl in `ls $input_path/*MT*R1*.gz | sed 's|.*/||'`; do
echo "Processing: "$fl
name=${fl::-12}
echo $name
echo "Saving at: "$outdir"/"$name
mkdir "$outdir"/"$name"
STAR --genomeDir /mnt/silencer2/share/STAR_indices/mm10 --readFilesIn $input_path/${name}_R1.fastq.gz $input_path/${name}_R2.fastq.gz --runThreadN 6 --outFileNamePrefix $outdir/$name/$name"_" --readFilesCommand zcat --quantMode GeneCounts
done
