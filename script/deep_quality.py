import json
from scipy.signal import find_peaks
import statistics
import pandas as pd
import sys
import subprocess
import os
from collections import defaultdict

def parsing(x, y):
    df = {}
    if y == "fastqc_per_base_sequence_quality_plot": 
        values = [] 
        for sample in x:
            for entry in sample:
                for data in entry['data']:
                    values.append(data[1])
        average_v = sum(values) / len(values)
        df["Mean per base sequence quality"] = average_v
    if y == "fastqc_per_sequence_quality_scores_plot":
        values = [] 
        for sample in x:
            for entry in sample:
                for data in entry['data']:
                    values.append(data[0])
        average = sum(values) / len(values)
        df["Mean per sequence quality scores"] = average
    if y == "fastqc_sequence_counts_plot":        
        values = []
        for sample in x:
            for entry in sample:
                    average = sum(entry['data']) / len(entry['data'])
                    df[entry['name']] = average
    if y == "fastqc_per_sequence_gc_content_plot":
        numb_peaks = []
        for sample in x:
            for entry in sample:
                values = [] 
                for data in entry['data']:
                    values.append(data[1])
                
                peaks, _ = find_peaks(values, height=4)
                numb_peaks.append(len(peaks)) 
        df["Mean number of peaks on GC plot"] = round(statistics.mean(numb_peaks))
    if y == "fastqc_per_base_n_content_plot":
        interval = []
        for sample in x:
            for entry in sample:
                nc = []
                for data in entry['data']:
                    if data[1] >= 6: 
                        nc.append(data[0])
                interval.extend(nc)    
        interval = sorted(set(interval))
        intervals = []
        if interval:
            start = interval[0]
            prev = interval[0]
            for val in interval[1:]:
                if val != prev + 2:
                    intervals.append([start, prev])
                    start = val
                prev = val
            intervals.append([start, prev])
        result = " ".join(f"[{x}-{y}]" if x != y else f"[{x}]" for x, y in intervals)
        df["Intervals of bases with N content"] = (result)
    if y == "fastqc_sequence_duplication_levels_plot":
        percent = []
        for sample in x:
            for entry in sample:
                percent.append(sum(entry['data'][9:]))
        df["Mean percent of libraries which have >10 duplications"] = statistics.mean(percent)
    if y == "fastqc_overrepresented_sequences_plot":
        values = []
        for sample in x:
            for entry in sample:
                    average = sum(entry['data']) / len(entry['data'])
                    df[entry['name']] = average
    else:
        pass
    
    return df

def generate_fastp_params(metrics_df, strictness='medium'):
    """
    Generates parameters for fastp based on MultiQC metrics.

    Input:
    metrics_df (pd.DataFrame): DataFrame with quality metrics
    strictness (str): Filter mode ('low', 'medium', 'high')

    Output:
    dict: Dictionary with parameters for fastp
    """
    metrics = dict(zip(metrics_df['Metrics'], metrics_df['Value']))
    # Default params
    params = {
        '--qualified_quality_phred': 20,
        '--unqualified_percent_limit': 40,
        '--n_base_limit': 5,
        '--dedup': '',
        '--correction': '',
        '--trim_front1': 0,
        '--trim_tail1': 0
    }

    # Setting the severity
    strictness_factor = {
        'low': 0.8,
        'medium': 1.0,
        'high': 1.2
    }[strictness]
    
    # Per base sequence quality
    mean_quality = metrics.get('Mean per base sequence quality', 30)
    params['--qualified_quality_phred'] = max(20, int(mean_quality * 0.7))
    params['--unqualified_percent_limit'] = int(50 / strictness_factor)
    
    # Duplicate
    if 'Unique Reads' in metrics and 'Duplicate Reads' in metrics:
        total_reads = metrics['Unique Reads'] + metrics['Duplicate Reads']
        if metrics['Duplicate Reads']/total_reads > 0.7:
            params['--dedup'] = '3'  # accuracy level to calculate duplication (1~6)
        else:
            params['--dedup'] = '1'
        
    # N-content
    if 'Intervals of bases with N content' in metrics:
        n_content = len(metrics['Intervals of bases with N content'].split())
        params['--n_base_limit'] = int(n_content * 0.5 * strictness_factor)
      
    # Overrepresented sequence
    if metrics.get('Top overrepresented sequence', 0) > 0.5:
        params['--overrepresentation_analysis'] = '--overrepresentation_analysis'  # enable overrepresented sequence analysis
        params['--overrepresentation_sampling'] = 20
        
    # Cut sequence 
    if mean_quality < 30:
        params['--trim_front1'] = 5
        params['--trim_tail1'] = 5
    
    return {k: v for k, v in params.items() if v != ''}

def process_fastp(args, pr):
    folder_path, filename = args
    if filename.endswith((".fastq", ".fq", "fastq.gz", "fq.gz")):
        input = folder_path + '/' + filename
        output = 'quality_after_filter/' + filename
        filter_seq = 'fastp -i {input} -o {output} {filt_params}'.\
            format(input=input, output=output, filt_params=' '.join('{0}'.format(w) for w in pr))
        exit_code = subprocess.call(filter_seq, shell=True)

def process_fastp_pair(args1, args2, folder, pr):
    if filename.endswith((".fastq", ".fq", "fastq.gz", "fq.gz")):
        input_1 = folder + '/' + args1
        input_2 = folder + '/' + args2
        output_1 = 'quality_after_filter/' + args1
        output_2 = 'quality_after_filter/' + args2
        print(input_1, input_2, output_1, output_2)
        filter_seq = 'fastp -i {input_1} -I {input_2} -o {output_1} -O {output_2} {filt_params}'.\
            format(input_1=input_1, input_2=input_2, output_1=output_1, output_2=output_2, filt_params=' '.join('{0}'.format(w) for w in pr))
        exit_code = subprocess.call(filter_seq, shell=True)

   
path = sys.argv[1]
fasta = sys.argv[2]

command = 'multiqc {dir} --outdir {dir} --interactive'.\
            format(dir=path)
exit_code = subprocess.call(command, shell=True)

with open(path + '/multiqc_data/multiqc_data.json', 'r') as file:
    data = json.load(file)

dataset = data['report_plot_data']['fastqc_per_base_sequence_quality_plot']['datasets']

all_mean = []
for w in data['report_plot_data'].keys():
    try:
        mean = parsing(data['report_plot_data'][w]['datasets'], w)
        for k, v in mean.items():
            all_mean.append({'Metrics': k, 'Value': v})
    except KeyError:
        pass

df_mean = pd.DataFrame(all_mean)    
if sys.argv[3] == 'on':
    params = json.loads(sys.argv[4])
    file_args = [(fasta, fn) for fn in os.listdir(fasta)]
    if sys.argv[6] == 'single':
       for w in file_args:
           results = process_fastp(w, params)
    if sys.argv[6] == 'paired':
       paired = defaultdict(dict)
       for w in file_args:
           folder_path, filename = w
           base = filename.replace("_1", "").replace("_2", "")
           if "_1" in filename:
               paired[base]["1"] = filename
           elif "_2" in filename:
               paired[base]["2"] = filename
       pairs = [(v["1"], v["2"]) for v in paired.values() if "1" in v and "2" in v]
       for r1, r2 in pairs:
           results = process_fastp_pair(r1, r2, folder_path, params)

if sys.argv[3] == 'off':
    fastp_params = generate_fastp_params(df_mean, strictness=sys.argv[5])

    print("Recommended parameters")
    params = []
    for k, v in fastp_params.items():
        print(f"{k} = {v}")
        params.append(f"{k} {v}")

    file_args = [(fasta, fn) for fn in os.listdir(fasta)]
#with mp.Pool(processes=2) as pool:
#        results = pool.map(partial(process_fastp, pr=params), file_args)
    if sys.argv[6] == 'single':
       for w in file_args:
           results = process_fastp(w, params)
    if sys.argv[6] == 'paired':
       paired = defaultdict(dict)
       for w in file_args:
           folder_path, filename = w
           base = filename.replace("_1", "").replace("_2", "")
           if "_1" in filename:
               paired[base]["1"] = filename
           elif "_2" in filename:
               paired[base]["2"] = filename
       pairs = [(v["1"], v["2"]) for v in paired.values() if "1" in v and "2" in v]
       for r1, r2 in pairs:
           results = process_fastp_pair(r1, r2, folder_path, params)



done = "filter done"
file2write=open("quality_after_filter/filter.txt",'w')
file2write.write(done)
file2write.close()
