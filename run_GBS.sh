#!/bin/bash

#SBATCH -A pronozinau
#SBATCH --mem 100GB
#SBATCH -p brain1
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -t 10:00:00

#python GBS_ID.py t.csv gunzip_fastq/
#source $CONDA_PREFIX/etc/profile.d/conda.sh
#conda activate GBS

snakemake -j 1
#Rscript cluster.R /hpcws/pronozinau-GBS/pipline/tree/cluster.gds