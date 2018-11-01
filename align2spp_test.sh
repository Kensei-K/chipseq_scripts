#!/bin/bash

PARENT_DIR=$PWD

# log files
LOG_FILE=$PARENT_DIR"/run_log.txt";LOG_ERR_FILE=$PARENT_DIR"/run_err.txt"
BARCODE_FILE=$PARENT_DIR/barcode.txt
SAMPLE_NO=`wc -l < $BARCODE_FILE`
NPROC_PER_SAMPLE=2
echo $SAMPLE_NO

cd 02trim

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.1: alignment" | tee -a $LOG_FILE

WORKING_DIR=$PARENT_DIR'/03alignment/'; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP

alignfun (){
    while read data;do
    samplename=${data%.*}
    echo -e "(`date`) align $samplename"
    bowtie2 -t --non-deterministic --mm --phred33 -p $NPROC_PER_SAMPLE --very-sensitive -x mm10 -U $data | samtools sort -@ $NPROC_PER_SAMPLE -T $samplename -o $WORKING_DIR"/"$samplename".bam" &
done
}
echo -e "(`date`) Performing alignment on files in here $PWD" | tee -a $LOG_FILE
ls -1 *.fastq | alignfun |tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP
echo -e "This version of bowtie2 was used for the analysis" | tee -a $LOG_FILE_STEP

bowtie2 —-version | tee -a $LOG_FILE_STEP 

echo -e "This version of samtools was used for the analysis" | tee -a $LOG_FILE_STEP

samtools —-version | tee -a $LOG_FILE_STEP 
wait;echo -e "`date`: Step 3.1 alignment Finished!" | tee -a $LOG_FILE


#------------------------------------------------------------
cd $WORKING_DIR
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.2: index bam file" | tee -a $LOG_FILE

date | tee -a $LOG_FILE_STEP
ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i\
                    samtools index {}\
                    1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP
wait;echo -e "(`date`) Step 3.2 index bam file Finished!" | tee -a $LOG_FILE


#------------------------------------------------------------
STEP='qualimap'
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR;

qualmapfun (){
    int1=1;int2=1;
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	mkdir -p $WORKING_DIR'/'$nameStr

	if [ `echo $int1" % 4" | bc` -eq 0 ]
	then
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G
	else
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G &
	fi
    done
}

ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
echo -e "This version of qualimap was used for the analysis" | tee -a $LOG_FILE_STEP
qualimap --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.1: filtering the aligned bam files" | tee -a $LOG_FILE
STEP='04filter'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i \
                     filterfun.sh {} $WORKING_DIR $NPROC_PER_SAMPLE \
    | tee -a $LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  Step 4.1 filter bam finished" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE

cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --progress -j $SAMPLE_NO samtools index {}
                    1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE

prefix=".txt"
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done | tee -a $LOG_FILE_STEP

date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  finished Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE




echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.3: Calculating NRF, PCR bottleneck coefficient" | tee -a $LOG_FILE


date | tee -a $LOG_FILE_STEP

# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

NRFqc () {
  while read data 
  do
	samplename=${data%.*}
	echo -e "(`date`) calculating NRF, PCR bottleneck coefficients for $samplename" | tee -a $LOG_FILE_STEP
	bedtools bamtobed -i $data | awk `BEGIN{OFS="\t"}{print $1,$2,$3,$6}` | grep -v `chrM` | sort | uniq -c | awk `BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}` > ${samplename}.txt 2>>$LOG_ERR_FILE
  done
}


ls *.bam | NRFqc | tee -a 1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE

date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  finished Step 4.3: Calculating NRF, PCR bottleneck coefficient" | tee -a $LOG_FILE









echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 5: deduplicate" | tee -a $LOG_FILE

STEP='05deDupped'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
mkdir -p $WORKING_DIR"/tmp"
#JAVA_OPTS=-Xmx64g
dedupfun (){
    while read data; do
    namestr=$(echo "$data"|cut -f1 -d".")
    java -Xmx64g -jar /opt/picard/picard.jar MarkDuplicates I=$data O=$WORKING_DIR"/"$namestr".bam" REMOVE_DUPLICATES=true ASSUME_SORTED=true M=$WORKING_DIR"/"$namestr".txt" TMP_DIR=$WORKING_DIR"/tmp" QUIET=F
done
}
ls -1 *.bam | dedupfun 1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
rm -r $WORKING_DIR"/tmp"
date | tee -a $LOG_FILE_STEP
echo -e "This version of MarkDuplicates from picard was used for the analysis" | tee -a $LOG_FILE_STEP
picard MarkDuplicates --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP


cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --progress -j $SAMPLE_NO samtools index {}
1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP

prefix=".txt"
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done | tee -a $LOG_FILE_STEP

date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  Step 4.3 deduplication finished" | tee -a $LOG_FILE






echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 6: spp cross-corelation " | tee -a $LOG_FILE
STEP='06spp'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"
LOG_ERR_FILE=$WORKING_DIR"/err.txt"

sppfun (){
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	Rscript /opt/phantompeakqualtools/run_spp.R -c="$data" -savp=$WORKING_DIR"/"$nameStr".pdf" -p=$NPROC_PER_SAMPLE &
    done
}

ls -1 *.bam | sppfun 1>> $LOG_FILE_STEP 2>> $LOG_ERR_FILE



STEP='qulimap'
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR;

qualmapfun (){
    int1=1;int2=1;
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	mkdir -p $WORKING_DIR'/'$nameStr

	if [ `echo $int1" % 4" | bc` -eq 0 ]
	then
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G
	else
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G &
	fi
    done
}

ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
echo -e "This version of qualimap was used for the analysis" | tee -a $LOG_FILE_STEP
qualimap --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP





##As the last step of this pipeline, plotFingerprint from deeptools is run
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 6.1: deduplicate" | tee -a $LOG_FILE
plotFingerprint --bamfiles *.bam --numberOfProcessors 2 -o fingerprint.png --outQualityMetrics fingerprintQC.txt | tee -a $LOG_FILE_STEP

echo -e "############################################################"
echo -e "This script has finished running"
