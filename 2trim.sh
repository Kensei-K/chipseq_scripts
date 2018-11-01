#!/bin/bash

echo -e "############################################################"
echo -e "(`date`) Welcome to the ChIPseq Auto analysis ver 0.0"
echo -e "############################################################"

echo -e "This pipeline will take a fastq file and finish it with filtered bam files"

echo -e "(`date`) Initiating input parameters"
echo -e " "

#input 1: consolidated FASTQ Folder
PARENT_DIR=$PWD #  project dir or current directory

# log files
LOG_FILE=$PARENT_DIR"/run_log.txt";LOG_ERR_FILE=$PARENT_DIR"/run_err.txt"

echo -e "(`date`) Setting up log files \n" | tee -a $LOG_FILE
echo -e "(`date`) Project folder is $PARENT_DIR \n" | tee -a $LOG_FILE
echo -e "(`date`) Log files are $LOG_ERR_FILE and $LOG_FILE \n" | tee -a $LOG_FILE
echo -e "------------------------------" | tee -a $LOG_FILE


echo -e "(`date`) Setting up log files \n" | tee -a $LOG_FILE
echo -e "(`date`) Project folder is $PARENT_DIR \n" | tee -a $LOG_FILE
echo -e "(`date`) Log files are $LOG_ERR_FILE and $LOG_FILE \n" | tee -a $LOG_FILE
echo -e "------------------------------" | tee -a $LOG_FILE
echo -e "(`date`) NOW please upload the password.txt file \n" | tee -a $LOG_FILE
echo -e "(`date`) It should only contain 1 line with : separator \n" | tee -a $LOG_FILE
echo -e "(`date`) example: SxaQSEQsWA146L6:wK4rq3yX4Em7 \n" | tee -a $LOG_FILE

read -p "Press [Enter] key to when you finished construct the file"

# barcode files
echo -e "------------------------------" | tee -a $LOG_FILE
echo -e "(`date`) NOW please upload the barcode.txt file \n" | tee -a $LOG_FILE
echo -e "(`date`) line format: samplename<tab>barcode \n" | tee -a $LOG_FILE
echo -e "(`date`) example: \n" | tee -a $LOG_FILE
echo -e "(`date`) example: \n" | tee -a $LOG_FILE
echo -e "(`date`) Wt-0	ATTCCTT
Wt-TNF-05	ACTGATA
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | te\Wt-TNF-1	GAGTGGA
Wt-TNF-3	CGTACGT
\n" | tee -a $LOG_FILE

read -p"Press [Enter] key to when you finished construct the file"


BARCODE_FILE=$PARENT_DIR/barcode.txt
echo -e "(`date`) barcode file is $BARCODE_FILE \n" | tee -a $LOG_FILE

# number of samples
SAMPLE_NO=`wc -l < $BARCODE_FILE`
echo -e "There are $SAMPLE_NO samples in this experiment: \n" | tee -a $LOG_FILE

while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    barcode=$(echo -e "$line" | cut -f2 -d$'\t')
    echo "$prefix $barcode"
    done<$BARCODE_FILE | tee -a $LOG_FILE


#number  of processors per sample for fastqc
NPROC_PER_SAMPLE=4
echo -e " $NPROC_PER_SAMPLE processors per sample \n" | tee -a $LOG_FILE
TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE)) # calculate number of total processor for the user
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | tee -a $LOG_FILE

echo -e "Check the password file (password.txt)" | tee -a $LOG_FILE
PASSWD_INFO_INPUT=$(head -n 1 password.txt)
echo -e "Your password is : $PASSWD_INFO_INPUT"

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) start running the pipelines " | tee -a $LOG_FILE
echo -e "(`date`) 1.1  starting downloading" | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

STEP="00raw"; WORKING_DIR=$PARENT_DIR"/"$STEP
LOG_FILE_STEP=$WORKING_DIR"/log.txt"
mkdir -p $WORKING_DIR; cd $WORKING_DIR
date| tee -a $LOG_FILE_STEP
grab_bscrc.sh $PASSWD_INFO_INPUT | tee - a $LOG_FILE_STEP

wait;echo -e "(`date`) downloaded the raw qseq data" | tee -a $LOG_FILE_STEP
RAW_DIR=$(echo $PASSWD_INFO_INPUT|cut -f1 -d":"); RAW_DIR=$WORKING_DIR/$RAW_DIR; cd $RAW_DIR
echo -e "(`date`) there are total: `ls -1 | wc -l`  raw qseq data" | tee -a $LOG_FILE_STEP

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.2 starting decompressing " | tee -a $LOG_FILE
date| tee -a $LOG_FILE_STEP
ls -1 *.gz | parallel -j $TOTAL_PROC_NO gunzip |tee -a $LOG_FILE
wait;echo -e "(`date`) decompressed the raw qseq data" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.3 starting converting qseq to fastq " | tee -a $LOG_FILE

##The first 10 letters of the qseq file will be used to make fastq files. 
qseq2fastqPar (){
    count=1;ncore=$TOTAL_PROC_NO;
    prefix=$WORKING_DIR"/"
    while read data
    do

        if [ `echo $count" % "$ncore | bc` -eq 0 ]
        then
	    qseq2fastq.pl $data ${data:0:10}'.fastq'
        else
	    qseq2fastq.pl $data ${data:0:10}'.fastq' &
        fi
        count=$((count+1))
    done
}

date| tee -a $LOG_FILE_STEP
ls -1 * | qseq2fastqPar | tee -a $LOG_FILE_STEP
ls -l | grep *.fastq | tee -a $LOG_FILE_STEP
wait;echo -e "(`date`) 1.3 converting finished" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e " (`date`) 1.4 starting demultiplex " | tee -a $LOG_FILE


STEP="consolidate";
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

demultiplexFunMuticore (){
    count=1
    ncore=$TOTAL_PROC_NO;
    prefix=$WORKING_DIR"/"
    while read data
    do
        suffix=${data:5:11}
        echo 'processing:'$suffix
        idxfile=$(echo "$data" |sed -r 's/_1_/_2_/')
        echo "using idx file: "$idxfile
        echo "calulating no.$count sampling ......"
	echo "suffix is $suffix"
	echo "save into folder $prefix"
	echo "(`date`) $count"
        if [ `echo $count" % "$ncore | bc` -eq 0 ] #mod ncore
        then
	    cat $data | fastx_barcode_splitter.pl --bcfile $BARCODE_FILE --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1
       else
           cat $data | fastx_barcode_splitter.pl --bcfile $BARCODE_FILE --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1 &
        fi
        count=$((count+1))
    done
}

date | tee -a $LOG_FILE_STEP
ls -1 *_1_*.fastq | demultiplexFunMuticore 1|tee -a $LOG_FILE_STEP 2>>$LOG_ERR_FILE
date | tee -a $LOG_FILE_STEP

wait;echo -e "(`date`) 1.4 demultiplex finished" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.5 start merge fastq " | tee -a $LOG_FILE


cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    eval "cat ${prefix}_* > ${prefix}.fastq"
    echo "(`date`)$prefix"
done<$BARCODE_FILE | tee -a $LOG_FILE_STEP 2>>$LOG_ERR_FILE
date | tee -a $LOG_FILE_STEP

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.6  start remove  original fastq " | tee -a $LOG_FILE
mkdir -p $WORKING_DIR'/unassigned'
mv unmatched* $WORKING_DIR'/unassigned'
rm *_*.fastq | tee -a $LOG_FILE_STEP
ls -lght >> $LOG_FILE




echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Part 1.7: qc of consolidated FASTQ files" | tee -a $LOG_FILE

WORKING_DIR=$WORKING_DIR'/fastqc';mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
		      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {}
                      1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
date | tee -a $LOG_ERR_FILE
echo -e "This version of fastqc was used for the analysis" | tee -a $LOG_FILE_STEP
fastqc --version >> $LOG_FILE_STEP 2>> $LOG_FILE_STEP
wait;echo -e "(`date`) Step 1.7 QC consolidated fastq Finished!" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2: trimming" | tee -a $LOG_FILE
WORKING_DIR=$PARENT_DIR'/02trim'
mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP  
ls -1 *.fastq | xargs -n1 -P $TOTAL_PROC_NO -i \
                      cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  {} \
                      -o $WORKING_DIR"/"{}".trim.fastq" \
                      1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
date | tee -a $LOG_FILE_STEP
echo -e "This version of cutadapt was used in the pipeline"
cutadapt --version | tee -a $LOG_FILE_STEP
echo -e "(`date`) Cut Adaptor of this sequence with, error rate of 0.1, phred quality score cut off was 20 and read shorter than 35 was discarded. Read and adapter have to overlap more than 6 bp to be trimmed." | tee -a $LOG_FILE
wait;echo -e "(`date`) Step 2: trimming  Finshed!"| tee -a $LOG_FILE



echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2.1: QC of trimed QC files" | tee -a $LOG_FILE
cd $WORKING_DIR
WORKING_DIR=$WORKING_DIR'/fastqc';mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i\
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {}\
                      1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
date | tee -a $LOG_ERR_FILE
wait;echo -e "(`date`)Step 2.1 QC trimed fastq Finshed!" | tee -a $LOG_FILE

