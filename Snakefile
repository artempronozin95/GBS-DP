configfile: "config.yaml"
import glob
import subprocess
import pandas as pd

dir = config['zip_fastq'].rsplit('/', 1)
if os.path.isdir(dir[0]) == True:
	files = glob.glob(config['zip_fastq'])
	zip = []
	gunzip = []
	for w in files:
    		out = w.rsplit('/', 1)
    		zip_out = out[1].rsplit('.', 1)
    		zip.append(out[1])
    		gunzip.append(zip_out[0])
	print(gunzip)
else:
	files = config['list_files']
	zip = []
	gunzip = []
	for w in files:
	        print(w)
    		out = w.rsplit('/', 1)
    		zip_out = out[1].rsplit('.', 1)
    		zip.append(out[1])
    		gunzip.append(zip_out[0])
	print(gunzip)

rule all:
    input:
        expand("{ref}.fa.bwt", ref=config['reference_genome']),
        expand("triming_reads/{sample}_trim.fastq", sample=gunzip),
        expand("alignment/{sample}_trim.fastq.sam", sample=gunzip),
        "alignment/res.txt",
        "alignment/sorted.txt",
        "VCF/vcf.txt",
        "mod.txt",
        "results/plotdendogram.png"

rule cutadupt:
    input:
        files
    output:
        expand("triming_reads/{sample}_trim.fastq", sample=gunzip)
    run:
        for w in input:
            print(w)
            out = w.rsplit('/', 1)
            out = out[1].rsplit('.', 1)
            shell("cutadapt -a AGATCGGAAGAG {w}  > triming_reads/{out[0]}_trim.fastq")

rule bwt:
   input:
           ref=ancient(expand("{ref}.fa", ref=config['reference_genome']))
   output:
           expand("{ref}.fa.bwt", ref=config['reference_genome'])
   run:
           shell("mkdir alignment")
           shell("bwa index {input.ref}")

rule alignment:
    input:
         trim_fastq=expand("triming_reads/{sample}_trim.fastq", sample=gunzip),
         index=expand("{ref}.fa.bwt", ref=config['reference_genome']),
         ref=expand("{ref}.fa", ref=config['reference_genome'])
    params:
            thr=4
    output:
           expand("alignment/{sample}_trim.fastq.sam", sample=gunzip)
    run:
        for w in input.trim_fastq:
            out = w.split('/')
            print(out)
            shell("bwa mem -t {params} {input.ref} {w} > alignment/{out[1]}.sam")

rule SAM_to_BAM:
    input:
            expand("alignment/{sample}_trim.fastq.sam", sample=gunzip)
    params:
           "alignment"
    output:
           "alignment/res.txt"
    run:
            shell("./script/sam_bam.sh {params}")
#            shell("rm -R triming_reads/")
#            shell("rm -R zip/")

rule sort:
    input:
        "alignment/res.txt"
    params:
        "alignment"
    output:
        "alignment/sorted.txt"
    run:
        shell("./script/sort.sh {params}")

rule VCF:
    input:
       "alignment/sorted.txt",
       ref=expand("{ref}.fa", ref=config['reference_genome'])
    params:
       fold1="alignment",
       fold2="VCF"
    output:
       "VCF/vcf.txt"
    run:
       shell("mkdir VCF_index")
       shell("script/VCF.sh {params.fold1} {params.fold2} {input.ref}")

rule index_vcf:
    input:
         vcf="VCF/vcf.txt",
         chr=expand("{chr}", chr=config['chr']),
         ref=expand("{ref}.fa", ref=config['reference_genome'])
    params:
         fold1="chr",
         fold2="VCF",
         fold3="VCF_index",
         fold4="tree"
    output:
         "mod.txt"
    run:
         shell("script/index_vcf.sh {params.fold1} {params.fold2} {params.fold3} {params.fold4} {input.ref}")

rule phyml:
    input:
         mod="mod.txt"
    params:
    	 vcf=ancient("tree/Merged.vcf"),
         fold1="chr/full_chr"
    output:
         "results/plotdendogram.png"
    run:
         shell("mkdir rtmp; TMPDIR=rtmp Rscript cluster_new.R {params.fold1} results/cluster.gds {input.mod} {params.vcf}")
