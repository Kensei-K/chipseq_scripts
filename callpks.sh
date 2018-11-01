#!/bin/bash


FRAG_LEN=200
LOG_FILE="run_log.txt";LOG_ERR_FILE="run_err.txt"

#Calling peaks for each sample against Input
############################################################
echo -e "(`date`) Starting Part 9: Calling peaks" | tee -a ../$LOG_FILE

mkdir 09MACS2
cd 05deDupped

callpeaks () {
  while read data 
  do
	samplename=${data%.*}
	echo -e "(`date`) calling peaks for ${samplename}"
	macs2 callpeak -t $data -c Input.bam -n $samplename -g mm -q 0.01 --keep-dup 1 --call-summits --nomodel \
	--extsize $FRAG_LEN --outdir ../09MACS2 >> ../09MACS2/$LOG_FILE 2>> ../09MACS2/$LOG_ERR_FILE
  done
}
date | tee -a ../09MACS2/$LOG_FILE
macs2 --version | tee -a ../09MACS2/$LOG_FILE
ls *.bam | grep -v 'I' | callpeaks |tee -a ../09MACS2/$LOG_FILE 2>> ../09MACS2/$LOG_ERR_FILE
