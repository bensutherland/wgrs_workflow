#!/bin/bash

# Global variables
GENOMEFOLDER="03_genome"
GENOME="GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz"
DATAFOLDER="04_samples"
NCPU="$1"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

# Align reads against the genome
for file in $(ls -1 "$DATAFOLDER"/*_R1.fastq.gz)
do
    # Name of second read file
    file2=$(echo "$file" | perl -pe 's/_R1.fastq.gz/_R2.fastq.gz/')
    echo "Aligning file $file $file2" 

    name=$(basename "$file")
    name2=$(basename "$file2")
    ID="@RG\tID:$name\tSM:$name\tPL:Illumina"

    # Align reads
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" "$DATAFOLDER"/"$name2" 2> /dev/null | 
        samtools view -Sb -q 10 - > "$DATAFOLDER"/"${name%.fastq.gz}".bam

    # Sort and index
    samtools sort --threads "$NCPU" -o "$DATAFOLDER"/"${name%.fastq.gz}".sorted.bam \
        "$DATAFOLDER"/"${name%.fastq.gz}".bam

    samtools index "$DATAFOLDER"/"${name%.fastq.gz}".sorted.bam

    # Cleanup
    rm "$DATAFOLDER"/"${name%.fastq.gz}".bam
done
