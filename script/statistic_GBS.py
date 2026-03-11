import pandas as pd
import glob
import subprocess
import re
import sys


gn_b = pd.read_csv(sys.argv[1], sep='\t')
gn_a = pd.read_csv(sys.argv[2], sep='\t')

gn_b = gn_b[['Sample','FastQC_mqc-generalstats-fastqc-percent_gc',
        'FastQC_mqc-generalstats-fastqc-avg_sequence_length',
        'FastQC_mqc-generalstats-fastqc-percent_duplicates']]

gn_a = gn_a[['Sample', 'FastQC_mqc-generalstats-fastqc-percent_duplicates']]

quality = pd.merge(gn_b, gn_a, how="inner", on=["Sample"])
quality.rename(columns={"FastQC_mqc-generalstats-fastqc-percent_gc" : "GC", 
 'FastQC_mqc-generalstats-fastqc-avg_sequence_length': "Sequence length",
 'FastQC_mqc-generalstats-fastqc-percent_duplicates_x': "Percent duplicates before", 
 'FastQC_mqc-generalstats-fastqc-percent_duplicates_y': "Percent duplicates after"}, inplace=True)

BAM = glob.glob("alignment/*.sort")

samples = []
reads = []
mapped_reads = []
mapped = []
depth = []
cover = []
for bam_file in BAM:

    cmd = f"samtools flagstat {bam_file}"
    result = subprocess.run(cmd, shell=True,  capture_output=True, text=True)
    output = result.stdout

    total_reads = int(re.search(r'^(\d+) \+ \d+ in total', output, re.MULTILINE).group(1))
    map_reads = int(re.search(r'(\d+) \+ \d+ mapped', output).group(1))
    mapping_rate = float(re.search(r'mapped \((\d+\.\d+)%', output).group(1))
    sample_id = re.search(r'/([^/]+?)(?:\..*)?$', bam_file).group(1)
    
    samples.append(sample_id)
    reads.append(total_reads)
    mapped_reads.append(map_reads)
    mapped.append(mapping_rate)

    bedtools = f"""bedtools genomecov -ibam {bam_file} -max 1 | awk '$1 == "genome" && $2 == 1 {{print $5}}'"""
    run_cover = subprocess.run(bedtools, shell=True, capture_output=True, text=True)
    output_cover = run_cover.stdout
    clean_output_cover = output_cover.strip()
    cover.append(clean_output_cover)
    
    samtools = f"samtools depth {bam_file} | awk '{{sum+=$3}} END {{print sum/NR}}'"
    result_depth = subprocess.run(samtools, shell=True, capture_output=True, text=True)
    mean_depth = float(result_depth.stdout.strip())
    depth.append(mean_depth)

BAM_data = pd.DataFrame(
    {'Sample': samples,
     'Number of reads': reads,
     'Number of mapped reads': mapped_reads,
     'Mapped, %':mapped,
     'Depth': depth,
     'Coverage': cover
    })

quality = pd.merge(quality, BAM_data, how="inner", on=["Sample"])
quality.to_csv(sys.argv[3], sep='\t', index=False)
