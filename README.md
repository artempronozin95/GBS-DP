# GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing

## Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Input](#input)
* [Running the pipeline](#work-start)
* [Output](#output)

## Introduction
The GBS method has demonstrated its reliability and flexibility for a number of plant species and populations. GBS has been applied to identify molecular markers for genetic mapping, genomic selection, genetic diversity research, variety identification, as well as research in conservation biology and evolutionary ecology. GBS method has reduced both the cost and the time required to sequencing of the studied samples. This has led to the need to develop high-quality bioinformatics analysis for an ever-expanding amount of sequenced data. For these purposes, to date, bioinformatics pipelines for analyzing data obtained by the GBS method have been developed.
In the precent work we developed bioinformatic pipeline GBS-DB. The pipeline is applicable for any species of organisms. The pipeline allows processing large amounts of data (more than 400 samples) and is implemented using Snakemake software manager

#### This pipeline is only applicable to the Linux operating system.

## Scheme of the bioinformatic pipeline
The pipeline includes the following steps: 
#### 1. Data pre-processing
+ Raw reads quaolity control
+ Adapters filtering
+ Reference genome index construction.
#### 2. Polymorphism searching
+ Reads mapping
+ Sorting of mapped reads
+ Polymorphism searching.
#### 3. Genetic diversity analysis
+ VCF format converting (into GDS)
+ Clasterization and phylogenetic tree construction.

The pipeline is implemented using the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), which provides ability to platform-independent installation and execution of the software.

## Schematic diagram
![Test Image 1](https://github.com/artempronozin95/GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing/blob/main/img/GBS_pipeline.png)

## Installation 
### Automatic
**recommended for clusters/servers**

Install only **programs.yaml** environment. Other environments will install automatically when ICAnnoLncRNA start work.
```
1. wget https://github.com/artempronozin95/GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing.git
2. unzip main.zip
3. cd ./GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing
4. conda env create --file env/programs.yaml
5. conda activate GBS
```
After these steps all necessary packages are installed. If you need update packages (**not recommended**), change the version of  packages after “=” (example - `snakemake=4.0.1 -> snakemake=6.0.0`), then `conda env update --file ./programs.yaml`. All necessary packages will be updated. Recomended on clusters, requires a lot of  processing power.
### Step method
**recommended for personal computer**
```
1. conda update conda.
2. conda create -n GBS python=3.6
3. conda activate GBS
4. conda install -c bioconda samtools
5. install next packeges from file below
```
## Input
### Raw sequence
Raw sequences obtened after sequencing methods in FASTQ.gz format.
### Reference genome
Reference genome of the species in `FASTA` format.

## Configuration file
Input all necessary files into configuration file “config.yaml”:
+ `zip_fastq:` - the path to the folder with the raw readings.
  + (Example: `zip_fastq: "zip/*.fastq.gz"`)
+ `reference_genome:` - the path to the reference genome and the name of genome.
  + (Example: `"ref/Prunus_persica_chr_numb"`)
    
## Work start
  #### 1. `snakemake -j 2`
  `-j` or  `--cores` -  Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
  #### 2. `snakemake -nr` 
  `-n` - Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.
  `-r` - Print the reason for each executed rule.
  #### 3. `snakemake --use-conda`
  `--use-conda` - Use additional conda environment.
  #### 4. Recommended run: 
  `snakemake -j 2 --use-conda`

## Output
A typical structure of `Output` is follows:
```
├── alignment
│   └── seq1.sam
│   └── seq1.bam
│   └── seq1.sort
│   └── seq2.sam
│   └── seq2.bam
│   └── seq2.sort
│   └── ...
├── chr
│   └── chr1
│   └── chr2
│   └── chr3
│   └── ...
├── env
│   └── programs.yaml
├── ref
│   └── reference.fasta
├── results
│   └── cluster.gds
│   └── cluster.png
│   └── cluster_tab.tsv
│   └── dendrogram.tree
│   └── plotdendogram.png
├── script
│   ├── index_vcf.sh
│   ├── sam_bam.sh
│   ├── sort.sh
│   └── VCF.sh
├── tree
│   └── Merged.vcf
│   └── Merged.vcf.stat
├── triming_reads
│   └── seq1.fastq
│   └── seq2.fastq
│   └── ...
├── VCF
│   └── seq1.vcf
│   └── seq2.vcf
│   └── ...
├── VCF_index
│   └── seq1.vcf
│   └── seq2.vcf
│   └── ...
└── zip
│   └── seq1.fastq.gz
│   └── seq2.fastq.gz
│   └── ...
```
***Data pre-processing***
+ zip - Raw reads in `FASTQ.gz` format.
+ triming_reads - processed reads, `FASTQ` format.
+ ref - reference genome, `FASTA` format. As well reference genome index in `FASTA.idx` format.
   
***Polymorphism searching***
+ alignment - reads mapping in `SAM`, `BAM` format. As well sorting mapping reads `SORT` format.
+ VCF - polymorphism searching, `VCF` format.
+ VCF_index -  index of polymorphism files in `CSI` format.
  
***Genetic diversity analysis***
+ chr - polymorphism file for each chromosome in `VCF` format. As well index of polymorphism files in `CSI` format, for each chromosome. 
+ tree - merged polymorphism files in `VCF` files.
+ results - containe main results of pipeline
  + cluster.gds - converted merged polymorphism files in `GDS` files.
  + cluster.png - PCA clusterization ![Test Image 2](https://github.com/artempronozin95/GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing/blob/main/img/cluster.png)
  + cluster_tab.tsv - Table of the PCA clusterization in `TSV` format.
  + dendrogram.tree - phylogenetic tree build by hierarchical clustering method in `TREE` format.
  + plotdendogram.png - phylogenetic tree build by hierarchical clustering method. ![Test Image 3](https://github.com/artempronozin95/GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing/blob/main/img/Tree.png)


