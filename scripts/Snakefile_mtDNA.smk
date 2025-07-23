configfile: "config_mtDNA.yaml"

#Use paths from config.yaml
READS_DIR = config["READS_DIR"]
OUTPUT_DIR = config["OUTPUT_DIR"]
REF_FASTA = config["REF_FASTA"]

import os
import random
from collections import defaultdict

# Get sample names based on folder structure in READS_DIR
SAMPLES = [d for d in os.listdir(READS_DIR) if os.path.isdir(os.path.join(READS_DIR, d))]

def get_paired_fastq_files(sample):
    """
    Efficiently match and randomly select paired-end FASTQ files.
    
    Args:
        sample (str): Sample name/directory
    
    Returns:
        tuple: Lists of matched R1 and R2 files with full paths
    """
    # Initialize dictionary to store file information
    file_pairs = defaultdict(dict)
    sample_path = os.path.join(READS_DIR, sample)
    
    # First pass: categorize files by their base name
    for filename in os.listdir(sample_path):
        if not filename.endswith('.fastq.gz'):
            continue
            
        if '_1.fastq.gz' in filename:
            base_name = filename.replace('_1.fastq.gz', '')
            file_pairs[base_name]['R1'] = os.path.join(sample_path, filename)
        elif '_2.fastq.gz' in filename:
            base_name = filename.replace('_2.fastq.gz', '')
            file_pairs[base_name]['R2'] = os.path.join(sample_path, filename)
    
    # Find complete pairs
    complete_pairs = [
        (pair['R1'], pair['R2'])
        for pair in file_pairs.values()
        if 'R1' in pair and 'R2' in pair
    ]
    
    # Randomly select pairs if we have more than 3
    if complete_pairs and len(complete_pairs) > 3:
        selected_pairs = random.sample(complete_pairs, 3)
    else:
        selected_pairs = complete_pairs
    
    # Split into R1 and R2 lists
    r1_files = [pair[0] for pair in selected_pairs]
    r2_files = [pair[1] for pair in selected_pairs]
    
    # Validation
    if not r1_files or not r2_files:
        raise ValueError(f"No complete pairs found for sample {sample}")
    
    return r1_files, r2_files

def get_fastq_files(wildcards):
    """Wrapper function for Snakemake rules"""
    r1_list, r2_list = get_paired_fastq_files(wildcards.sample)
    return {'R1': r1_list, 'R2': r2_list}

# Rule to define the final outputs of the workflow
rule all:
    input:
        expand(f"{OUTPUT_DIR}/consensus/{{sample}}_mtDNA_consensus.fasta", sample=SAMPLES)

rule index_reference:
    input:
        ref = f"{OUTPUT_DIR}/{REF_FASTA}"
    output:
        indexes = expand(f"{OUTPUT_DIR}/{REF_FASTA}.{{ext1}}", ext1=["bwt.2bit.64", "pac", "0123", "ann", "amb"]),
        fai = f"{OUTPUT_DIR}/{REF_FASTA}.fai",
    shell:
        """
        export PATH=/programs/bwa-mem2-2.2.1:$PATH
        bwa-mem2 index {input.ref}
        samtools faidx {input.ref} --output {output.fai}
        """

rule merge_fastq:
    input:
        unpack(get_fastq_files)
    output:
        temp_R1=temp(f"{OUTPUT_DIR}/temp/{{sample}}_temp_R1.fastq"),
        temp_R2=temp(f"{OUTPUT_DIR}/temp/{{sample}}_temp_R2.fastq")
    shell:
        """
        # Concatenate R1 and R2 separately and unzip
        zcat {input.R1:q} > {output.temp_R1}
        zcat {input.R2:q} > {output.temp_R2}
        """

rule fastq_pair:
    input:
        temp_R1=f"{OUTPUT_DIR}/temp/{{sample}}_temp_R1.fastq",
        temp_R2=f"{OUTPUT_DIR}/temp/{{sample}}_temp_R2.fastq"
    output:
        pair_R1=temp(f"{OUTPUT_DIR}/temp/{{sample}}_temp_R1.fastq.paired.fq"),
        pair_R2=temp(f"{OUTPUT_DIR}/temp/{{sample}}_temp_R2.fastq.paired.fq")
    shell:
        """
        # Use fastq-pair to synchronize paired-end reads
        /programs/fastq_pair/bin/fastq_pair -t 3109179 {input.temp_R1} {input.temp_R2}
        """

rule fastq_to_fastq_gz:
    input:
        pair_R1=f"{OUTPUT_DIR}/temp/{{sample}}_temp_R1.fastq.paired.fq",
        pair_R2=f"{OUTPUT_DIR}/temp/{{sample}}_temp_R2.fastq.paired.fq"
    output:
        R1=temp(f"{OUTPUT_DIR}/merged/{{sample}}_merged_1.fastq.gz"),
        R2=temp(f"{OUTPUT_DIR}/merged/{{sample}}_merged_2.fastq.gz")
    threads: 32  # Number of threads to use for pigz
    shell:
        """
        # pigz fastq
        /programs/pigz-2.4/pigz -p {threads} -c {input.pair_R1} > {output.R1}
        /programs/pigz-2.4/pigz -p {threads} -c {input.pair_R2} > {output.R2}
        """

rule fastq_to_bam:
    input:
        ref = REF_FASTA,
        R1=f"{OUTPUT_DIR}/merged/{{sample}}_merged_1.fastq.gz",
        R2=f"{OUTPUT_DIR}/merged/{{sample}}_merged_2.fastq.gz",
        indexes = expand(f"{OUTPUT_DIR}/{REF_FASTA}.{{ext}}", ext=["bwt.2bit.64", "pac", "0123", "ann", "amb", "fai"]),
    output:
        bam = f"{OUTPUT_DIR}/bams/{{sample}}.bam",
    log:
        f"{OUTPUT_DIR}/bwa_mem_log/{{sample}}.txt"
    shell:
        """
        export PATH=/programs/bwa-mem2-2.2.1:$PATH
        bwa-mem2 mem -t 64 {input.ref} {input.R1} {input.R2} 2> {log} | samtools sort -T {wildcards.sample} -o {output.bam}

        """

rule index_bam:
    input:
        bam = f"{OUTPUT_DIR}/bams/{{sample}}.bam"
    output:
        bai = f"{OUTPUT_DIR}/bams/{{sample}}.bam.bai",
    shell:
        """
        samtools index {input.bam} > {output.bai}
        """

rule bam_to_vcf:
    input:
        bam = f"{OUTPUT_DIR}/bams/{{sample}}.bam",
        ref = REF_FASTA,
    output:
        vcf = f"{OUTPUT_DIR}/vcfs/{{sample}}_variants.vcf.gz"
    shell:
        """
        bcftools mpileup -Oz -f {input.ref} {input.bam} | \
        bcftools call --ploidy 1 -mv -Oz -o {output.vcf} && \
        bcftools index -t {output.vcf}
        """

rule consensus_mtDNA:
    input:
        ref = REF_FASTA,
        vcf = f"{OUTPUT_DIR}/vcfs/{{sample}}_variants.vcf.gz"
    output:
        output = f"{OUTPUT_DIR}/consensus/{{sample}}_mtDNA_consensus.fasta"
    shell:
        """
        bcftools consensus -f {input.ref} {input.vcf} > {output.output}
        sed -i 's/^>/>{wildcards.sample}_/' {output.output}
        """

        