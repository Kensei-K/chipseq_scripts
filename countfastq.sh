#!/bin/bash

echo "Provide a relative path to where your fastq files are"
read dir

for filename in $dir/*.fastq; do
	reads=$(awk '{s++}END{print s/4}' $filename)
	echo "$filename has $reads reads"
done

