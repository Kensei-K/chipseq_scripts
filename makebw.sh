#!/bin/bash

echo "Provide a relative path to where your bam files"
read dir

for filename in $dir/*.bam; do
	name=$(echo $filename | sed " s@$dir@@")
	name=$(echo $name | sed " s@/@@")
	name=$(echo $name | sed " s@.bam@@")
	echo "Start processing $name" | tee -a ./BWLOG.txt
        bamCoverage --bam $filename -o ./$name.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX --extendReads 200 --numberOfProcessors 6 -v | tee -a ./BWLOG.txt
	echo "Done processing $name"
done
