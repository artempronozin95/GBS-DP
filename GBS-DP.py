#!/usr/bin/env python3

# Python Standard Packages
import argparse
import subprocess
import os
import sys

# Other Packages
import pandas as pd
from yaml import dump as yaml_dump


# %% Functions

def create_config(configfile, args):
    #endregion # settings
    #region # create config dictionary
    args.snakefile == os.path.join(os.path.dirname(os.path.realpath(__file__)), "Snakefile")
    ref_id = args.reference_genome.rsplit('.', 1)
    print(args.sample_files)
    config_dict = {
        'list_files' : args.sample_files,
        'zip_fastq'   : os.path.join(args.sample_dir, '*.fastq.gz'),
        'reference_genome'    : ref_id[0],
        'chr' :   args.chromosomes  }

    # endregion # create config dictionary

    # write dictionary to yaml
    with open(configfile, 'w') as outfile:
        yaml_dump(config_dict, outfile, default_flow_style = False, sort_keys = False)
    print(f"Config written to {configfile}")


def run_snakemake(configfile, args):

    # other workflow options
    threads = "--cores" if args.threads == "none" else ""  # use all available cores, if no profile is specified
    threads = f"--cores {args.threads}" if args.threads > 0 else threads  # override cores by commandline arg if provided
    dryrun = "--dryrun --quiet" if args.dryrun else ""
    unlock = "--unlock" if args.unlock else ""

    call = f"snakemake --configfile {configfile} " \
           f"--snakefile {args.snakefile}  {threads} {dryrun} {unlock}"

    if args.use_conda:
        frontend = "conda" if args.conda_frontend else "mamba"
        call = call + f" --use-conda --conda-prefix {args.condaprefix} --conda-frontend {frontend}"

    print(call)
    subprocess.call(call, shell = True)

# %% Main

def main():
    # set global vars
    global pipeline
    global repo_path
    pipeline = "GBSDP"

    # set paths
    repo_path = os.path.dirname(os.path.realpath(__file__))

    #region # parse arguments
    parser = argparse.ArgumentParser()

    # path arguments
    parser.add_argument('--sample_files', '-l', help = '[REQUIRED] list of samples to analyze, fastq files.',
                        default = False,  type=str, nargs='+', required = False)
    parser.add_argument('--sample_dir', '-d', help = '[REQUIRED] directory of samples to analyze, fastq file paths.',
                        default = "/", type = os.path.abspath, required = False)
    parser.add_argument('--reference_genome', '-r', help = '[REQUIRED] Reference genome of the studied organism',
                        default = False, type = os.path.abspath, required = False)
    parser.add_argument('--chromosomes', '-chr', help = '[REQUIRED] name of chromosomes of the studied organism',
                        default = False, type = os.path.abspath, required = False)
    parser.add_argument('--snakefile', '-s', help = f'Path to Snakefile of {pipeline}; default: {os.path.join(repo_path, "Snakefile")}',
                        default = os.path.join(repo_path, "Snakefile"), type = os.path.abspath, required = False),
    # cpu arguments
    parser.add_argument('--use_conda', help = 'Utilize the Snakemake "--useconda" option, i.e. Smk rules require execution with a specific conda env',
                        default = False, action = 'store_true', required = False)
    parser.add_argument('--conda_frontend', help = 'Do not use mamba but conda as frontend to create individual conda environments',
                        default = False, action = 'store_true', required = False)
    parser.add_argument('--condaprefix', '-c', help = f'Path of default conda environment, enables recycling of built environments; default: {os.path.join(repo_path, "conda_env")}',
                        default = os.path.join(repo_path, "conda_env"), required = False)
    parser.add_argument('--threads', '-t', help = f'Number of Threads/Cores to use. This overrides the "<{pipeline}>/profiles" settings',
                        default = 0, type = int, required = False)
    parser.add_argument('--dryrun', '-n', help = 'Snakemake dryrun. Only calculate graph without executing anything',
                        default = False, action = 'store_true', required = False)
    parser.add_argument('--unlock', help = 'Unlock a Snakemake execution folder if it had been interrupted',
                        default = False, action = 'store_true', required = False)



    args = parser.parse_args()

       # checks
#    if args.sample_dir == "/":
#        print("ERROR: the following argument is required: -d/--sample_dir")
#        sys.exit(1)
    if not os.path.exists(args.sample_dir):
        print(f"ERROR: path to sample list {args.sample_dir} does not exist")
        sys.exit(1)

    if args.reference_genome == False:
        print("ERROR: the following argument is required: -r/--reference_genome")
        sys.exit(1)


    #endregion # eval arguments

    # write config
    configfile = "./config.yaml"
    create_config(configfile, args)

    run_snakemake(configfile, args)




# %% Main Call

if __name__ == '__main__':
    main()
