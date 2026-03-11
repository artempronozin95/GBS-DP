import subprocess
import tempfile
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# HIST
def hist(path):
    # Read tables
    afreq = pd.read_csv(f'{path}/st_merged_bisnp.afreq', sep='\t')
    vmiss = pd.read_csv(f'{path}/st_merged_bisnp.vmiss', sep='\t')
    smiss = pd.read_csv(f'{path}/st_merged_bisnp.smiss', sep='\t')
    hwe = pd.read_csv(f'{path}/st_merged_bisnp.hardy', sep='\t')
    log = open(f'{path}/LOG_stats.txt', "w")

    # Add total call count
    vmiss['TOTAL'] = vmiss['MISSING_CT'] + vmiss['OBS_CT']
    vmiss['F_MISS_CALC'] = vmiss['MISSING_CT'] / vmiss['TOTAL']

    # ALT allele frequency histogram
    plt.figure(figsize=(8,6))
    sns.histplot(afreq['ALT_FREQS'], bins=50, color='skyblue', edgecolor='black')
    plt.title("Distribution of ALT Allele Frequencies")
    plt.xlabel("ALT Allele Frequency")
    plt.ylabel("Count of Variants")
    plt.tight_layout()
    plt.savefig(f'{path}/bisnp_altfreq.png')
    plt.close()


    # Minor allele frequency
    afreq['MAF'] = np.minimum(afreq['ALT_FREQS'], 1 - afreq['ALT_FREQS'])
    plt.figure(figsize=(8,5))
    sns.histplot(afreq['MAF'], bins=50, color='skyblue', edgecolor='black', alpha=0.8)
    plt.axvline(0.01, color='red', linestyle='--', lw=1)
    plt.axvline(0.05, color='orange', linestyle='--', lw=1)
    plt.text(0.015, plt.gca().get_ylim()[1]*0.9, "MAF = 0.01", color="red", fontsize=8)
    plt.text(0.065, plt.gca().get_ylim()[1]*0.9, "MAF = 0.05", color="orange", fontsize=8)
    plt.title("MAF Distribution with Common Filters")
    plt.xlabel("Minor Allele Frequency (MAF)")
    plt.ylabel("Count")
    plt.xlim([0,0.5])
    plt.xticks(np.arange(0, 0.51, 0.05))
    plt.tight_layout()
    plt.savefig(f'{path}/bisnp_maf.png', dpi=300)
    plt.close()


    # Per-variant missingness
    plt.figure(figsize=(8,6))
    sns.histplot(vmiss['F_MISS'], bins=50, color='steelblue', edgecolor='black', alpha=0.8)
    plt.title("Distribution of Per-Variant Missingness (F_MISS)")
    plt.xlabel("Missingness Rate")
    plt.ylabel("Number of Variants")
    plt.xlim([0, vmiss['F_MISS'].max()])
    plt.tight_layout()
    plt.savefig(f'{path}/bisnp_vmiss.png')
    plt.close()


    # Per-individual missingness
    plt.figure(figsize=(8,6))
    sns.histplot(smiss['F_MISS'], bins=30, color='green', edgecolor='black')
    plt.title("Individual Genotyping Missing Rate")
    plt.xlabel("Fraction Missing (per individual)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f'{path}/bisnp_smiss.png')
    plt.close()


    # Hardy-Weinberg histogram
    hwe['P'] = hwe['P'].replace(0, 1e-300)
    log_p_values = np.log10(hwe['P'])
    plt.figure(figsize=(8,6))
    sns.histplot(log_p_values, bins=50, color='purple', edgecolor='black')
    min_log = np.floor(log_p_values.min())
    max_log = np.ceil(log_p_values.max())
    ticks = np.arange(min_log, max_log + 1, step=10)
    plt.xticks(
    ticks=ticks,
    labels=[f'$10^{{{int(x)}}}$' for x in ticks])
    plt.title("HWE P-values (Unaffected Controls)")
    plt.xlabel("HWE P-value (log scale)")
    plt.ylabel("Variant Count")
    plt.tight_layout()
    plt.savefig(f'{path}/bisnp_hwe.png')
    plt.close()

    # Summary statistics
    print("Summary statistics of MAF", file=log)
    print(afreq['MAF'].describe(), file=log)
    print("Number of variants with MAF < 0.01:", (afreq['MAF'] < 0.01).sum(), file=log)
    print("Proportion of variants with MAF > 0.05 (common variants):", (afreq['MAF'] > 0.05).mean(), file=log)
    print("Summary of missingness per variant", file=log)
    print(vmiss['F_MISS'].describe(), file=log)
    print("Mean missing rate per variant:", vmiss['F_MISS'].mean(), file=log)
    print("Median missing rate per variant:", vmiss['F_MISS'].median(), file=log)
    print("Min missing rate:", vmiss['F_MISS'].min(), file=log)
    print("Max missing rate:", vmiss['F_MISS'].max(), file=log)
    high_miss = (vmiss['F_MISS'] > 0.2).sum()
    print(f"Variants with >20% missingness: {high_miss} ({100*high_miss/len(vmiss):.2f}%)", file=log)

# PCA
def PCA_plot(list, path):
    for w in list:
        pca = pd.read_csv(f'{path}/{w}.eigenvec', sep='\t')
        eigenval = np.loadtxt(f'{path}/{w}.eigenval')

        pca = pca.rename(columns={'#IID': 'ind'})
     
        pve = pd.DataFrame({
            'PC': np.arange(1, 11),
            'pve': (eigenval[:10] / eigenval.sum()) * 100})
        
        plt.figure(figsize=(7,5))
        sns.barplot(x='PC', y='pve', data=pve, color='skyblue', edgecolor='black')
        plt.ylabel('Percentage variance explained')
        plt.xlabel('PC')
        plt.tight_layout()
        plt.savefig(f'{path}/{w}_barplot.png')
        plt.close()

        plt.figure(figsize=(7,7))
        sns.scatterplot(data=pca, x='PC1', y='PC2', s=90, edgecolor=None, palette='Paired')
        plt.xlabel(f"PC1 ({pve['pve'][0]:.3g}%)")
        plt.ylabel(f"PC2 ({pve['pve'][1]:.3g}%)")
        plt.gca().set_aspect('equal', adjustable='datalim')
        sns.despine()
        plt.grid(False)
        plt.tight_layout()
        plt.savefig(f'{path}/{w}_PCA.png')
        plt.close()

def main():
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    missing = sys.argv[3]
    maf = sys.argv[4]
    mind = sys.argv[5]
    LD_act = sys.argv[6]

    merged_true_path = f'{output_path}/filtered.vcf'

    # Фильтрация строк и запись непосредственно в файл в output_path
    with open(input_path, 'r') as f_in, open(merged_true_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                parts = line.strip().split('\t')
                if len(parts) >= 5 and parts[3] != '.' and parts[4] != '.':
                    f_out.write(line)

# biallelic
    
    biallelic = 'plink2 --vcf {merged_true_path} --make-pfile --chr-set 50 --allow-extra-chr --set-all-var-ids @:# \
                 --max-alleles 2 --snps-only just-acgt --min-alleles 2 --out {output_path}/merged_bisnp' .\
            format(merged_true_path=merged_true_path, output_path=output_path)
    
    exit_biallelic = subprocess.call(biallelic, shell=True)
    

    biallelic_2 = 'plink2 --pfile {output_path}/merged_bisnp  --freq --missing --hardy --out {output_path}/st_merged_bisnp' .\
        format(output_path=output_path)
    
    exit_biallelic_2 = subprocess.call(biallelic_2, shell=True)
    
# Hist build
    hist(output_path)

    TOM_hardy = pd.read_csv(f'{output_path}/st_merged_bisnp.hardy', sep='\t')
    total = TOM_hardy['HOM_A1_CT'] + TOM_hardy['HET_A1_CT'] + TOM_hardy['TWO_AX_CT']
    TOM_hardy['het_rate'] = TOM_hardy['O(HET_A1)']/total
    
    TOM_hardy = TOM_hardy[TOM_hardy['het_rate'] < 0.10]
    
    TOM_hardy['ID'].to_csv(f'{output_path}/keep_variants.txt', sep='\t', index=False)
# missing data
 
    missing = 'plink2 --pfile {output_path}/merged_bisnp --chr-set 50 --allow-extra-chr --set-all-var-ids @:# \
                --extract {output_path}/keep_variants.txt --geno {missing} --export vcf --out {output_path}/bisnp_het.miss2' .\
            format(output_path=output_path, missing=missing)
    exit_missing = subprocess.call(missing, shell=True)
   # os.remove(merged_true_path)


# Individs
    Individs = 'plink2 --vcf {output_path}/bisnp_het.miss2.vcf --chr-set 50 --allow-extra-chr --set-all-var-ids @:# \
                       --mind {mind} --export vcf --out {output_path}/bisnp_het.miss2.ind' .\
            format(output_path=output_path, mind=mind)
    exit_missing = subprocess.call(Individs, shell=True)


# MAF

    MAF = 'plink2 --vcf {output_path}/bisnp_het.miss2.ind.vcf --chr-set 50 --allow-extra-chr --set-all-var-ids @:# \
                  --maf {maf} --export vcf --out {output_path}/bisnp_het.miss2.maf01' .\
            format(output_path=output_path, maf=maf)
    exit_missing = subprocess.call(MAF, shell=True)


# LD
    if LD_act == 'on':
        r = sys.argv[7]
        LD = 'plink2 --vcf {output_path}/bisnp_het.miss2.maf01.vcf --chr-set 50 --allow-extra-chr --set-all-var-ids @:# \
                       --rm-dup force-first --indep-pairwise 50 10 {r} --out {output_path}/bisnp_het.miss2.maf01.ind.LD' .\
            format(output_path=output_path, r=r)
        exit_missing = subprocess.call(LD, shell=True)
        
        SNP = 'plink2 --vcf {output_path}/bisnp_het.miss2.maf01.vcf --chr-set 50 --allow-extra-chr --set-all-var-ids @:# \
                       --extract {output_path}/bisnp_het.miss2.maf01.ind.LD.prune.in --export vcf --out {output_path}/{output_prefix}' .\
            format(output_path=output_path, output_prefix='bisnp_het.miss2.maf01.ld_on')
        exit_missing = subprocess.call(SNP, shell=True)
    else:
        pass
# PCA table
    if LD_act == 'on':
        list = ['bisnp_het.miss2', 'bisnp_het.miss2.maf01', 'bisnp_het.miss2.maf01.ld_on']
    else:
    #    list = ['bisnp_het.miss2', 'bisnp_het.miss2.maf01']
        list = ['bisnp_het.miss2.maf01']
    for w in list:
        PCA = 'plink2 --vcf {output_path}/{file}.vcf --pca 12 --out {output_path}/{file}' .\
            format(output_path=output_path, file=w)
        exit_missing = subprocess.call(PCA, shell=True)
    PCA_plot(list, output_path)
    
if __name__ == '__main__':
    main()
