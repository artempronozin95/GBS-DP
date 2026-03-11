cat $1/vcf.txt | parallel -j 40 --tmpdir ./ "vcffilter -f 'MQ > $3' -f 'DP > $4' -f 'QUAL > $5' $1/{} > $2/{}"
#| awk '\$5 != \"\\.\"'> $2/{}"
ls $2/*.vcf | awk -F '/' '{print $2}' > $2/vcf.txt
