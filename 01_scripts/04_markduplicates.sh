#!/bin/bash
# Mark duplicates in each BAM (PCR or optical)

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Set variables
JAVA="/usr/bin/java" # full path needed
PICARD="/usr/local/bin/picard-tools-1.131/picard.jar"        # full path needed

# Mark duplicates on all sorted bam files 
for file in $(ls -1 04_samples/*.sorted.bam)
do
    # Reporting
    echo "Marking duplicates for: " $file

    # Calculate coverage
    $JAVA -jar $PICARD MarkDuplicates I=$file \
        O="${file%.bam}"_mdups.bam \
        M="${file%.bam}"_mdups_metrics.txt 
        #--OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 

done 2>&1 | tee 10_log_files/"$TIMESTAMP"_picard_mdups.log



