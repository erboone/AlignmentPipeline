# TODO: Look into the best way to read .gz files for alignment
# TODO: add robustness of fastq/fasta, etc in input
# TODO: Any reason to keep the samfiles around
    # Look into overhead of bam > sam


rule all:
    input:
        "data/bigwig/{sample}.bigwig"

rule _directory_setup_model:
    output:
        '.chkpts/mkdir_chkpt'
    params:
        directories=directory(expand("data/{subdir}", 
                              subdir=["raw", "aligned", "logs", "bigWig"]))
    shell:
        "mkdir -p {params.directories}"
        "touch .chkpts/.mkdir_chkpt"

# TODO: come back to the downloading part later
"""rule download_sra_target:
    input:
        expand('data/raw/{sample}_{pair_no}.fastq' ,
                pair_no = ["1","2"])"""

rule _download_sra_model:
    output:
        "data/raw/{sample}_1.fastq",
        "data/raw/{sample}_2.fastq"
    shell:
        """
        echo "Downloading accession ID: {wildcards.sample}"
        fasterq-dump --split-files {wildcards.sample} -O data/raw
        gzip {output}
        """

rule StarAlign:
    input:
        ".chkpts/mkdir_chkpt",
        ref_genome="data/ref_genome.fa",
        seq_data1="data/raw/{sample}_1.fastq.gz",
        seq_data2="data/raw/{sample}_2.fastq.gz"

    output:
        aligned="data/aligned/{sample}.sam"
        #TODO: Add QC output

    # TODO: look into making this customizeable, or at least raise it for better computers
    threads:2
    shell:
        """
        echo "Processing: {wildcards.sample}"
        echo "Saving at: {output.aligned}"
        STAR --genomeDir {input.ref_genome} \
             --readFilesIn {input.seq_data1} {input.seq_data2}\
             --runThreadN {threads} \
             --outFileNamePrefix {output.aligned} \
             --readFilesCommand zcat \
             --quantMode GeneCounts
        """
    

rule SamToBam:
    input:
        "data/aligned/{sample}.sam"
    output:
        "data/aligned/{sample}.sorted.bam"
    shell:
        """
        samtools view -bS {input} | samtools sort -o {output}
        
        """
        # TODO: remove SAM? Combine this with the StarAlignRule?

rule BamToBigWig:
    input:
        "data/aligned/{sample}.sorted.bam"
    output:
        "data/bigwig/{sample}.bigwig"
    shell:
        """
        echo "Converting: "{wildcards.sample}" to bigWig file."
        echo "Saving as:'{output}'"
        bedtools genomeCoverageBed -d -bga -g ./ -ibam {input}  > {output}
        """