#!/bin/bash

echo "Provide a relative path to where your bam files are"
read dir

for filename in $dir/*.bam; do
        reads=$(samtools view -c -F 260 $filename)
        echo "$filename has $reads reads"
done

