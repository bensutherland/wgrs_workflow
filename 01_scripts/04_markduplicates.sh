#!/bin/bash
# Mark duplicates in each BAM (PCR or optical)

# Set variables
JAVA="/home/greent/programs/jdk-21.0.1/bin/java" # full path needed
PICARD="/home/greent/programs/picard.jar"        # full path needed

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

done

