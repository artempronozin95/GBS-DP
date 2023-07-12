cat $1/sorted.txt| parallel -j 40 --tmpdir ./ "samtools mpileup -u -t AD -aa --fasta-ref $3 $1/{} | bcftools call -m -o $2/{}.vcf"
ls $2/*.vcf | awk -F '/' '{print $2}' > $2/vcf.txt