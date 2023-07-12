# GBS-DP-bioinformatics-pipeline-for-genotyping-by-sequencing-data-processing


## Introduction
#### This pipeline is only applicable to the Linux operating system.


The pipeline includes the following steps: 


The pipeline is implemented using the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), which provides ability to platform-independent installation and execution of the software.

## Schematic diagram


## Installation 
# Automatic
**recommended for clusters/servers**

Install only **programs.yaml** environment. Other environments will install automatically when ICAnnoLncRNA start work.
```
1. wget https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/archive/refs/heads/main.zip
2. unzip main.zip
3. cd ./ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA-main
4. conda env create --file env/programs.yaml
5. conda activate ICAnnoLncRNA
```
After these steps all necessary packages are installed. If you need update packages (**not recommended**), change the version of  packages after “=” (example - `snakemake=4.0.1 -> snakemake=6.0.0`), then `conda env update --file ./programs.yaml`. All necessary packages will be updated. Recomended on clusters, requires a lot of  processing power.
# Step method
**recommended for personal computer**
```
1. conda update conda.
2. conda create -n ICAnnoLncRNA python=3.6
3. conda activate ICAnnoLncRNA
4. conda install -c bioconda emboss
5. install next packeges from file below
```
## Input



### Configuration file

    
### Folders

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
├── results
│   └── cluster.gds
│   └── cluster.png
│   └── cluster_tab.tsv
│   └── dendrogram.tree
│   └── plotdendogram.png
├── rtmp
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
```

  

Load environment:
1. conda env create --file env/programs.yaml
2. conda activate GBS

Input files that should wright into config:
1. folder with .gz fastq files. (sample: zip/*.fastq.gz) 
2. reference genome.
All other job will do pipeline.
https://docs.google.com/document/d/1wGryeajmjzpX6ysa2g-yuacTAYoZzM_cIFIIvPGVKjw/edit#
