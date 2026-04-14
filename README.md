# Genotype a known deletion with amplicon sequencing

Jun Ishigohoka


## Description

This Snakemake pipeline genotypes deletions based on amplicon sequencing data.
The amplicon should be designed around the deletion breakpoint.

### Input

- Paired-end amplicon reads in FASTQ
- `{locus}.fa`: Reference of focal region in FASTA
- A space delimited info table 
    - column 1: sample
    - column 2: locus
    - column 3: allele
    - column 4: primer
- `{locus}_{allele}_del.bed`: A BED file specifying positions of deletion breakpoints
- `{locus}_{allele}_flank.bed`: A BED file specifying flanking regions of the deletion (to compute background depth)



### Output

- `trimmed/`: Directory containing trimmed FASTQ files
- `bam/{sample}/`: Directory containing read mapping.
- `genotype/genotype.txt`: A text file with six columns (sample, locus, allele, genotype, depth_del, depth_flank)


### Workflow

1. Trimming of reads using `trim_galore`
2. Index reference FASTA (for `bwa-mem2` later)
3. Map trimmed reads to reference using `bwa-mem2`, generating a BAM file
4. Extract soft-clipped reads from the BAM file using `samtools` and a custom `awk`
5. Re-map extracted soft-clipped reads to reference using `bwa-mem2`
6. Merge inital and remapped BAM files
7. Call genotype based on depth inside and outside the deletion breakpoints

### Dependencies

- `conda`
- `snakemake`

Other dependencies should be resolved by `conda` while running `snakemake`


## Run pipeline


### Prepare reference file in FASTA


Obtain the reference sequence in FASTA e.g. from NCBI.


### Prepare BED files


A BED file is a tab-delimited table with 3 columns specifying genomic regions.
- column 1: sequence name (e.g. chromosome)
- column 2: start position (0-based, inclusive)
- column 3: end position 2 (1-based, exclusive)

For example, chr_1 11-20 and chr_2 1-15 are written


```
chr_1	10	20
chr_2	0	15
```


You need to prepare two BED files, one specifying the delesion, and another specifying flanking sequences of the delesion.
The sequence name (1st column) should match the sequence name that appears in the reference FASTA, and the coordinate should be based on this sequence.
For example, if the deletion is chr_1:1000101-1000200, but the reference sequence starts from chr_1:1000001, then the 2nd and 3rd columns should be 100 and 200.
Also note that if you use gene sequence from NCBI, the orientation is opposite when the locus is on the complementary side of the reference genome.


The BED files have to be named in the following format:

- `{sample}_{locus}_{allele}_del.bed`
- `{sample}_{locus}_{allele}_flank.bed`

### Prepare info table

Prepare a space delimited table with 4 columns.
Currently, space and underscore (_) are not allowed.

- column 1: sample
- column 2: locus
- column 3: target allele
- column 4: primer


### Edit `config.yaml`

The input files/folders need to exist.
Output folders will be created based on the config.



## Run the pupeline

Also check <https://github.com/PallaresLab/DNAPipeline> and <https://snakemake.readthedocs.io/en/stable/> to learn how to run `snakemake`.

```bash
snakemake --snakefile Snakefile --configfile ./config.yaml --use-conda --conda-prefix ./envs/ --cores 1

```


