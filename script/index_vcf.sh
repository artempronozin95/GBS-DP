size=`du -c $2 | tail -1 | cut -f 1`
echo $size
if [[ $size -gt 100000000000 ]]
then
mkdir chr
echo 'the total volume of all VCF file is higher than 1 Tb. Use chunk сhromosome method'
cat $2/vcf.txt | parallel -j 20 --tmpdir ./ "bgzip $2/{} -c > $3/{.}.gz ; bcftools index $3/{.}.gz"

ls $3/*.gz | awk -F '/' '{print $2}' > $3/vcf.txt

for w in $(cat ref/chr.txt)
do
mkdir $1/$w
done

for w in $(cat $3/vcf.txt)
do
echo $w
cat ref/chr.txt| parallel -j 20 --tmpdir ./ "bcftools view -Oz $3/$w --regions {} -o $1/{.}/$w ; bcftools index $1/{.}/$w"
done

mkdir $1/full_chr

cat ref/chr.txt | parallel -j 20 --tmpdir ./ "bcftools merge -Ov $1/{}/*.gz -o $1/full_chr/{.}.vcf" 

ls $1/full_chr/*.vcf | awk -F '/' '{print $3}' > $1/full_chr/allchr.txt

cat $1/full_chr/allchr.txt | parallel -j 20 --tmpdir ./ 'cat chr/full_chr/{} | awk '\''$0 ~ /^#/ || ($4 != "." && $5 != ".") {print $0}'\'' > chr/full_chr/{}.filt'

mkdir tree

parallel bgzip {} ::: $1/full_chr/*.filt
parallel bcftools index ::: $1/full_chr/*.filt.gz

bcftools concat $1/full_chr/*.filt.gz -Ov -a -d all -o $4/Merged.vcf
bcftools stats -F $5 -s - $4/Merged.vcf > $4/Merged_stat.vcf

python script/vcf_filter.py $4/Merged.vcf $4 $6 $7 $8 $9 $10

echo 'chunk_mod' > mod.txt

else
echo 'the total volume of all VCF file is less than 1 Tb. Use simple method'
mkdir tree
cat $2/vcf.txt | parallel -j 20 --tmpdir ./ "bgzip $2/{} -c > $3/{.}.gz ; bcftools index $3/{.}.gz"
bcftools merge -Ov $3/*.gz -o $4/Merged.vcf
bcftools stats -F $5 -s - $4/Merged.vcf > $4/Merged_stat.vcf
ls $3/*.gz | awk -F '/' '{print $2}' > $3/vcf.txt

python script/vcf_filter.py $4/Merged.vcf $4 $6 $7 $8 $9 $10

echo 'simple_mod' > mod.txt
fi
