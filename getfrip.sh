#!/bin/sh


echo "Please enter the name of a txt file that names of bams files to be processes"
read bam

echo "Please enter the name of a bed file that contains peaks"
read bed

plotEnrichment -b $(less $bam) --BED $bed -v --smartLabels -e 250 -o frip.png -p 6 --variableScales --outRawCounts frip.tab| tee -a log.txt


