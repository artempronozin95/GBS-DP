configfile: "config.yaml"
import glob
import subprocess
import pandas as pd
from collections import defaultdict

files = glob.glob(config['zip_fastq'])
zip = []
gunzip = []
#if files.endswith((".fastq", ".fq", "fastq.gz", "fq.gz")):
for w in files:
    out = w.rsplit('/', 1)
    zip_out = out[1].split('.')
    zip.append(out[1])
    gunzip.append(zip_out[0])

rule all:
    input:
        expand("{ref}.bwt", ref=config['reference_genome']),
        expand("quality_after_filter/{sample}_fastqc.html", sample=gunzip),
        expand("quality_before_filter/{sample}_fastqc.html", sample=gunzip),
        "quality_after_filter/multiqc_report.html",
        "quality_after_filter/filter.txt",
        "alignment/res.txt",
        "alignment/sorted.txt",
        "VCF/vcf.txt",
        "VCF_filt/vcf.txt",
        "mod.txt",
        "results/plotdendogram.png",
        "results/lib_stat.csv"

rule fastqc:
   input:
       files
   params:
       type=expand("{type}", type=config['type'])
   output:
       expand("quality_before_filter/{sample}_fastqc.html", sample=gunzip)
   run:
    if params.type[0] == 'single':
       for w in input:
            print(w)
            shell("fastqc {w} -o quality_before_filter/")
    if params.type[0] == 'paired':
       paired = defaultdict(dict)
       for f in input:
            base = f.replace("_1", "").replace("_2", "")
            if "_1" in f:
               paired[base]["1"] = f
            elif "_2" in f:
               paired[base]["2"] = f
       pairs = [(v["1"], v["2"]) for v in paired.values() if "1" in v and "2" in v]
       for r1, r2 in pairs:
           shell("fastqc {r1} {r2} -o quality_before_filter/")


rule filter:
   input:
       fastqc=ancient(expand("quality_before_filter/{sample}_fastqc.html", sample=gunzip)),
       fasta=expand("{folder}", folder=config['folder'])
   params:
       type=expand("{type}", type=config['type']),
       multiqc="quality_before_filter",
       option=expand("{option}", option=config['fastp']['option']),
       metrics=lambda wildcards: json.dumps([f'{k} {v}' for k, v in config['fastp']['metrics'].items()]),
       strictness_factor = expand("{strictness_factor}", strictness_factor=config['fastp']['strictness_factor'])
   output:
       "quality_after_filter/filter.txt"
   conda:
       "env/filter.yaml"
   shell:
       "python script/deep_quality.py {params.multiqc} {input.fasta} {params.option} '{params.metrics}' {params.strictness_factor} {params.type}"

rule fastqc_after:
   input:
     "quality_after_filter/filter.txt"
   params:
     type=expand("{type}", type=config['type']),
     files="quality_after_filter/*.fq.gz"
   output:
     expand("quality_after_filter/{sample}_fastqc.html", sample=gunzip)
   run:
      files = glob.glob(str(params.files))
      if params.type[0] == 'single':
         for w in files:
            print(w)
            shell("fastqc {w} -o quality_after_filter/")
      if params.type[0] == 'paired':
         paired = defaultdict(dict)
         for f in files:
            base = f.replace("_1", "").replace("_2", "")
            if "_1" in f:
               paired[base]["1"] = f
            elif "_2" in f:
               paired[base]["2"] = f
         pairs = [(v["1"], v["2"]) for v in paired.values() if "1" in v and "2" in v]
         for r1, r2 in pairs:
            shell("fastqc {r1} {r2} -o quality_before_filter/")

rule multiqc:
   input:
     expand("quality_after_filter/{sample}_fastqc.html", sample=gunzip)
   output:
     "quality_after_filter/multiqc_report.html"
   conda:
     "env/filter.yaml"
   shell:
     "multiqc quality_after_filter/ --outdir quality_after_filter/"

rule bwt:
   input:
       ref=ancient(expand("{ref}", ref=config['reference_genome']))
   output:
       expand("{ref}.bwt", ref=config['reference_genome'])
   run:
       shell("mkdir alignment")
       shell("bwa index {input.ref}")

rule alignment:
  input:
       start=ancient("quality_after_filter/filter.txt"),
       index=expand("{ref}.bwt", ref=config['reference_genome']),
       ref=expand("{ref}", ref=config['reference_genome'])
  params:
       type=expand("{type}", type=config['type']),
       trim_fastq=expand("quality_after_filter/{sample}.fq.gz", sample=gunzip),
       thr=4
  output:
       "alignment/res.txt"
  run:
     if params.type[0] == 'single':
       for w in params.trim_fastq:
            out = w.split('/')
            print(out)
            shell("bwa mem -t {params.thr} {input.ref} {w} | samtools sort -o alignment/{out[1]}.bam -")
       shell("""ls alignment/*.bam | awk -F '/' '{{print $2}}' > alignment/res.txt""")
     if params.type[0] == 'paired':
       paired = defaultdict(dict)
       for f in params.trim_fastq:
           base = f.replace("_1", "").replace("_2", "")
           if "_1" in f:
               paired[base]["1"] = f
           elif "_2" in f:
               paired[base]["2"] = f
       pairs = [(v["1"], v["2"]) for v in paired.values() if "1" in v and "2" in v]
       for r1, r2 in pairs:
           out = r1.split('/')
           out = out[1].split('_')
           print(out)
           shell("bwa mem -t {params.thr} {input.ref} {r1} {r2} | samtools sort -o alignment/{out[0]}.bam -")
       shell("""ls alignment/*.bam | awk -F '/' '{{print $2}}' > alignment/res.txt""")

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
       ancient("alignment/sorted.txt"),
       ref=expand("{ref}", ref=config['reference_genome'])
   params:
       fold1="alignment",
       fold2="VCF"
   output:
       "VCF/vcf.txt"
   run:
#       shell("mkdir VCF_index")
#       shell("mkdir VCF_filt")
       shell("script/VCF.sh {params.fold1} {params.fold2} {input.ref}")

rule filter_VCF:
    input:
       vcf="VCF/vcf.txt"
    params:
       fold1="VCF",
       fold2="VCF_filt",
       MQ=expand("{MQ}", MQ=config['VCF_filter']['MQ']),
       QUAL=expand("{QUAL}", QUAL=config['VCF_filter']['QUAL']),
       DP=expand("{DP}", DP=config['VCF_filter']['DP'])
    output:
       "VCF_filt/vcf.txt"
    run:
       shell("script/VCF_filt.sh {params.fold1} {params.fold2} {params.MQ} {params.DP} {params.QUAL}")

rule idex_vcf:
    input:
         vcf="VCF_filt/vcf.txt",
         chr="ref/chr.txt",
         ref=expand("{ref}", ref=config['reference_genome'])
    params:
         fold1="chr",
         fold2="VCF_filt",
         fold3="VCF_index",
         fold4="tree",
         maf=expand("{maf}", maf=config['VCF_filter']['MAF']),
         missing=expand("{MS}", MS=config['VCF_filter']['Missing_data']),
         mind=expand("{MIND}", MIND=config['VCF_filter']['mind']),
         LD_act=expand("{LD_act}", LD_act=config['VCF_filter']['act']),
         LD=expand("{LD}", LD=config['VCF_filter']['LD'])
    output:
         "mod.txt"
    run:
         shell("script/index_vcf.sh {params.fold1} {params.fold2} {params.fold3} {params.fold4} {input.ref} {params.missing} {params.maf} {params.mind} {params.LD_act} {params.LD}")

rule phyml:
    input:
         mod="mod.txt"
    params:
         vcf=ancient("tree/bisnp_het.miss2.maf01.vcf"),
         fold1="chr/full_chr"
    output:
         "results/plotdendogram.png"
    run:
         shell("mkdir rtmp; TMPDIR=rtmp Rscript cluster_new.R {params.fold1} results/cluster.gds {input.mod} {params.vcf}")

rule statistics:
    input:
       "quality_after_filter/multiqc_report.html"
    params:
       after="quality_after_filter/multiqc_data/multiqc_general_stats.txt",
       before="quality_before_filter/multiqc_data/multiqc_general_stats.txt",
    output:
       "results/lib_stat.csv"
    conda:
       "env/stat.yaml"
    shell:
       "python script/statistic_GBS.py {params.before} {params.after} {output}"
