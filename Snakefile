# TODO: Look into the best way to read .gz files for alignment
# TODO: add robustness of fastq/fasta, etc in input
# TODO: Any reason to keep the samfiles around
    # Look into overhead of bam > sam

SAMPLES = ["SRR27488309"]
LOG_DIR = "data/_log/"
LOG_FILES = ["0_setup", "1QC_fastqc", "1_Hisat2Align", "1_SamToBam", "A2_BamToBigWig", "B2_FeatureCounts"]


rule all:
    input:
        expand("data/counts/{samples}_counts.txt",samples=SAMPLES),
        expand("data/_log/.1qc.{samples}.done",samples=SAMPLES),

# -----------------------------------------------------------------------------
# Setup:
# The first step of the pipeline, should be run before anything else: sets up 
# directory infreastructure
# TODO: Download necessary packages 
# TODO: create an updated DAG of jobs as an aid in documentation

rule:
    name: "0_setup"
    output:
        touch(f"{LOG_DIR}0_setup.log")
    log:
        f"{LOG_DIR}0_setup.log"
    params:
        directories=directory(expand("data/{subdir}", 
                              subdir=["_log", "raw", "reference", "1_aligned", "A2_bigwig", 'B2_counts'])),
    shell:
        """
        mkdir -p {params.directories} 
        """


# -----------------------------------------------------------------------------
# Alignment:
# The primary function of this snakemake pipeline. Aligns raw NGS data to ref. 
# genome. Additional quality control functionality.
# directory infreastructure

rule:
    name: "1QC_fastqc"
    input:
        seq_data1="data/raw/{sample}_1.fastq.gz",
        seq_data2="data/raw/{sample}_2.fastq.gz"
    output:
        touch(f"{LOG_DIR}.1qc."+ "{sample}.done"),
        fq1="data/1_aligned/{sample}_1_fastqc",
        fq2="data/1_aligned/{sample}_2_fastqc"

    script:
        """
        fastqc -o {output.fq1} {input.seq_data1}
        fastqc -o {output.fq2} {input.seq_data2}
        """

# Intermediate 
rule: 
    name: "1_Hisat2Align"
    input:
        ref_genome="data/reference/ref_genome.fa",
        seq_data1="data/raw/{sample}_1.fastq.gz",
        seq_data2="data/raw/{sample}_2.fastq.gz"

    output:
        aligned=temp("data/1_aligned/{sample}.sam")
    log:
        f"{LOG_DIR}" + "{sample}_1_Hisat2Align.log"
    threads:2
    shell:
        """
        echo "Aligning Reads from {wildcards.sample} using Hisat2"
        hisat2  -x {input.ref_genome} \
                -1 {input.seq_data1} \
                -2 {input.seq_data2} \
                --rna-strandness RF \
                -S {output.aligned} \
                --threads {threads}

        """

rule:
    name: "1_SamToBam"
    input:
        "data/1_aligned/{sample}.sam"
    output:
        "data/1_aligned/{sample}.sorted.bam"
    log:
        f"{LOG_DIR}" + "{sample}_1_SamToBam.log"

    shell:
        """
        samtools view -bS {input} | samtools sort -o {output}
        """
        # TODO: remove SAM? Combine this with the StarAlignRule?

rule: 
    name: "A2_BamToBigWig"
    input:
        "data/1_aligned/{sample}.sorted.bam"
    output:
        bw_plus ="data/A2_bigwig/{sample}_plus.bw",
        bw_minus="data/A2_bigwig/{sample}_minus.bw"
    log: 
        f"{LOG_DIR}" + "{sample}_A2_BamToBigWig.log"
    shell:
        """
        echo "Converting: "{wildcards.sample}" to bigWig file."
        echo "Saving as:'{output}'"
        bamCoverage --filterRNAstrand forward --normalizeUsing CPM --effectiveGenomeSize 2913022398  -b {input} -o {output.bw_plus}
        bamCoverage --filterRNAstrand reverse --normalizeUsing CPM --effectiveGenomeSize 2913022398  -b {input} -o {output.bw_minus}
        """
        #""" OLD IMPLEMENTATION: Consider deleting
        #echo "Converting: "{wildcards.sample}" to bigWig file."
        #echo "Saving as:'{output}'"
        #bedtools genomeCoverageBed -d -bga -g ./ -ibam {input}  > {output}
        #"""

rule: 
    name: "B2_FeatureCounts"
    input:
        bam="data/1_aligned/{sample}.sorted.bam",
        gtf="data/reference/features.gtf"
    output:
        "data/counts/{sample}_counts.txt"
    log: 
        f"{LOG_DIR}" + "{sample}_B2_FeatureCounts.log"
    script:
        """
        htseq-count -f bam \
                    -s reverse \
                    -t transcript \
                    -r pos \
                    -m intersection-nonempty \
                    --nonunique all \
                    -t transcript \
                    -i gene_name {input.bam} \
                     {input.gtf}
                    > {output}
        """


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# These are targets that I implemented earlier, but need to be polished to 
# productively fit into the worflow.

"""rule download_sra_target:
    input:
        expand('data/raw/{sample}_{pair_no}.fastq' ,
                pair_no = ["1","2"])

    rule _download_sra_model:
    output:
        "data/raw/{sample}_1.fastq",
        "data/raw/{sample}_2.fastq"
    shell:

        echo "Downloading accession ID: {wildcards.sample}"
        fasterq-dump --split-files {wildcards.sample} -O data/raw
        gzip {output}
"""


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
        echo "Processing: {wildcards.sample}"
        echo "Saving at: {output.aligned}"
        STAR --genomeDir {input.ref_genome} \
             --readFilesIn {input.seq_data1} {input.seq_data2}\
             --runThreadN {threads} \
             --outFileNamePrefix {output.aligned} \
             --readFilesCommand zcat \
             --quantMode GeneCounts
"""

