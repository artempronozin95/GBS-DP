cat $1/res.txt | parallel -j 40 --tmpdir ./ samtools sort $1/{} -o $1/{.}.sort
ls $1/*.sort | awk -F '/' '{print $2}' > $1/sorted.txt