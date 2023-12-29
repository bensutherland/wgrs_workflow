#!/bin/bash
# Remove adapters using cutadapt
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
NCPU=$1

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"
ADAPTERS_FILE="00_archive/adapters.fasta" 

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

# Create directory for untrimmed files
mkdir 02_raw/trimmed 2>/dev/null

rm 10_log_files/"$TIMESTAMP"_01_cutadapt"${i%.fastq.gz}".log 2> /dev/null

# Paired-end mode
ls -1 02_raw/*R1_001.fastq.gz | 
        sed 's/R1_001\.fastq\.gz//g' | 
        sed 's/02_raw\///g' | 
parallel -j $NCPU cutadapt -a file:"$ADAPTERS_FILE" \
        -A file:"$ADAPTERS_FILE" \
        -o 02_raw/trimmed/{}"R1.fastq.gz" \
        -p 02_raw/trimmed/{}"R2.fastq.gz" \
        -e 0.2 \
        -q 15 \
        --trim-n \
        -m 50 \
        "02_raw/"{}"R1_001.fastq.gz" "02_raw/"{}"R2_001.fastq.gz" '2>&1' '>>' 10_log_files/"$TIMESTAMP"_01_cutadapt.log


