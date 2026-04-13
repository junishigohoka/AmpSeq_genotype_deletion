import glob
import os
import pandas as pd

configfile: "config.yaml"


dir_fastq_raw = config["dir_fastq_raw"]
dir_analysis = config["dir_analysis"]
dir_fastq_trimmed = dir_analysis + "/trimmed"
dir_bam = dir_analysis + "/bam"
dir_genotype = dir_analysis + "/genotype"
dir_ref = config["dir_ref"]
dir_bed = config["dir_bed"]
n_threads = config["n_threads"]

samples_df = pd.read_table(config["info_table"], sep='\s+', header=None, names=["id", "locus", "allele", "primer"])

wildcard_constraints:
    sample = "[^/_]+",
    locus = "[^/_]+",
    allele = "[^/_]+"

def get_paired_fastq(wildcards):
        R1 = glob.glob(dir_fastq_raw + f"/{wildcards.sample}*_R1_001.fastq.gz")
        R2 = glob.glob(dir_fastq_raw + f"/{wildcards.sample}*_R2_001.fastq.gz")
        if len(R1) != 1:
                raise ValueError(f"R1 not found for {wildcards.sample}")
        if len(R2) != 1:
                raise ValueError(f"R2 not found for {wildcards.sample}")
        return [R1[0], R2[0]]





rule all:
        input:
                dir_genotype + "/genotypes.txt"

rule trim:
        input:
                get_paired_fastq
        output:
                R1 = dir_fastq_trimmed + "/{sample}/{sample}_{locus}_{allele}_val_1.fq.gz",
                R2 = dir_fastq_trimmed + "/{sample}/{sample}_{locus}_{allele}_val_2.fq.gz",
        log:
                dir_fastq_trimmed + "/{sample}/{sample}_{locus}_{allele}.log"
        params:
                dir_out = dir_fastq_trimmed + "/{sample}",
                basename = "{sample}_{locus}_{allele}"
        threads:
                n_threads
        conda:
                "envs/trim.yaml"
        shell:
                """
                mkdir -p {params.dir_out}
                trim_galore --paired {input}  --basename {params.basename} --fastqc --output_dir {params.dir_out} --cores {threads} >{log} 2>&1
                """



rule index_ref:
        input: 
                ref = dir_ref + "/{locus}.fa"
        output: 
                fa = temp(dir_bam + "/{locus}.fa"),
                bwt = temp(dir_bam + "/{locus}.fa.bwt.2bit.64"),
                idx = temp(dir_bam + "/{locus}.fa.0123"),
                amb = temp(dir_bam + "/{locus}.fa.amb"),
                ann = temp(dir_bam + "/{locus}.fa.ann"),
                pac = temp(dir_bam + "/{locus}.fa.pac")
        params:
                dir_bam = dir_bam
        log:
                dir_bam + "/{locus}_bwa_index.log"
        conda:
                "envs/bwa.yaml"
        shell:
                """
                mkdir -p {params.dir_bam}
                ln -s $(realpath {input}) {output.fa}
                bwa-mem2 index {output.fa} >{log} 2>&1
                """

rule map_reads:
        input:
                R1 = dir_fastq_trimmed + "/{sample}/{sample}_{locus}_{allele}_val_1.fq.gz",
                R2 = dir_fastq_trimmed + "/{sample}/{sample}_{locus}_{allele}_val_2.fq.gz",
                bwt = dir_bam + "/{locus}.fa.bwt.2bit.64",
                idx = dir_bam + "/{locus}.fa.0123",
                amb = dir_bam + "/{locus}.fa.amb",
                ann = dir_bam + "/{locus}.fa.ann",
                pac = dir_bam + "/{locus}.fa.pac",
                ref = dir_bam + "/{locus}.fa"
        output:
                bam = dir_bam + "/{sample}/{sample}_{locus}_{allele}.bam",
                bai = dir_bam + "/{sample}/{sample}_{locus}_{allele}.bam.bai"
        params:
                dir_out = dir_bam + "/{sample}"
        threads: 
                n_threads
        conda:
                "envs/bwa.yaml"
        shell:
                """
                mkdir -p {params.dir_out}
                bwa-mem2 mem -M -t {threads} {input.ref} {input.R1} {input.R2} | samtools sort -@{threads} | samtools view -b -F 256 -f 2 -o {output.bam}
                samtools index {output.bam}
                """

rule extract_clip:
        input: dir_bam + "/{sample}/{sample}_{locus}_{allele}.bam"
        output: dir_bam + "/{sample}/{sample}_{locus}_{allele}_softclip.fa"
        shell:
                """
                samtools view {input} | awk '$6 ~ /^[0-9]+S/' | awk '{{
                match($6, /^([0-9]+)S/, arr)
                print ">"$1"\\n"substr($10, 1, arr[1])
                }}' > {output}
                """

rule map_clipped_reads:
        input: 
                reads = dir_bam + "/{sample}/{sample}_{locus}_{allele}_softclip.fa",
                ref = dir_bam + "/{locus}.fa",
                bwt = dir_bam + "/{locus}.fa.bwt.2bit.64",
                idx = dir_bam + "/{locus}.fa.0123",
                amb = dir_bam + "/{locus}.fa.amb",
                ann = dir_bam + "/{locus}.fa.ann",
                pac = dir_bam + "/{locus}.fa.pac",
                bam_ori = dir_bam + "/{sample}/{sample}_{locus}_{allele}.bam"
        output: 
                bam_remapped = dir_bam + "/{sample}/{sample}_{locus}_{allele}_softclip_remapped.bam",
                bai_remapped = dir_bam + "/{sample}/{sample}_{locus}_{allele}_softclip_remapped.bam.bai",
                bam_combined = dir_bam + "/{sample}/{sample}_{locus}_{allele}_combined.bam",
                bai_combined = dir_bam + "/{sample}/{sample}_{locus}_{allele}_combined.bam.bai"
        conda:
                "envs/bwa.yaml"
        shell:
                """
                bwa-mem2 mem {input.ref} {input.reads} | samtools sort -o {output.bam_remapped}
                samtools index {output.bam_remapped}
                samtools merge -f {output.bam_combined} {input.bam_ori} {output.bam_remapped} 
                samtools index {output.bam_combined}
                """

rule call_genotype:
    input:
        bed_del = dir_bed + "/{locus}_{allele}_del.bed",
        bed_flank = dir_bed + "/{locus}_{allele}_flank.bed",
        bam_combined = dir_bam + "/{sample}/{sample}_{locus}_{allele}_combined.bam"
    output:
        temp(dir_genotype + "/{locus}/{allele}/{sample}_{locus}_{allele}_genotype.txt")
    params:
        dir_out = dir_genotype + "/{locus}/{allele}"
    conda:
        "envs/call_genotype.yaml"
    shell:
        """
        mkdir -p {params.dir_out}
        depth_del=`bedtools coverage -d -a {input.bed_del} -b {input.bam_combined} | awk '{{s+=$5;i++}}END{{print s/i}}'`
        depth_flank=`bedtools coverage -d -a {input.bed_flank} -b {input.bam_combined} | awk '{{s+=$5;i++}}END{{print s/i}}'`
        awk -v id={wildcards.sample} -v locus={wildcards.locus} -v allele={wildcards.allele} -v depth_del=$depth_del -v depth_flank=$depth_flank 'BEGIN{{if(depth_del < depth_flank / 4){{gt=2}}else if(depth_del > 3 * depth_flank / 4){{gt=0}}else{{gt=1}};printf "%s %s %s %d %.2f %.2f\\n", id, locus, allele, gt, depth_del, depth_flank}}' > {output}
        """

rule aggregate:
    input:
        expand(dir_genotype + "/{locus}/{allele}/{sample}_{locus}_{allele}_genotype.txt",
               zip,
               sample=samples_df["id"],
               locus=samples_df["locus"],
               allele=samples_df["allele"])
    output:
        dir_genotype + "/genotypes.txt"
    shell:
        "cat {input} > {output}"
