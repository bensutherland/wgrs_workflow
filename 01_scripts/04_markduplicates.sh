#!/bin/bash
# Mark or remove duplicates in each BAM (PCR or optical)
#  note: set REM_DUPS=false if you want to only mark duplicates

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Set variables
JAVA="/usr/bin/java" # full path needed
PICARD="/home/greent/programs/picard.jar"        # full path needed
REM_DUPS=true

# Mark duplicates on all sorted bam files 
for file in $(ls -1 04_samples/*.sorted.bam)
do
    # Reporting
    echo "Marking duplicates for: " $file

    # Mark duplicates, and remove (if specified by variable)
    $JAVA -jar $PICARD MarkDuplicates I=$file \
        O="${file%.bam}"_mdups.bam \
        M="${file%.bam}"_mdups_metrics.txt \
        REMOVE_DUPLICATES=$REM_DUPS

done 2>&1 | tee 10_log_files/"$TIMESTAMP"_picard_mdups.log

