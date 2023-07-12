ls $1/*.sam | parallel -j 40 --tmpdir ./ samtools view -S -b -o {.}.bam {}
ls $1/*.bam | awk -F '/' '{print $2}' > $1/res.txt
