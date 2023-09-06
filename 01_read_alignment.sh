#!/bin/bash

# Set up the running directory, raw data and reference genome.
export RawData="/home/boffelli/mystorage/Data"
export Running="/home/boffelli/mystorage/running"
export REF="$RawData/RefGenome/hg19.fa"

# Load the reference in the memory (only run this once) 
bwa shm $REF
# Run the alignment for parent, pacc, and progeny.
for file in $RawData/fastq/{controls,pacc,progeny}/*.fastq.gz; do
    dir="${file%/*}"
    dir="${dir##*/}"
    sample=${file##*/}
    bwa mem -R "@RG\tID:$sample\tSM:${dir/s/}\tLB:${dir/s/}\tPL:ILLUMINA\tPU:ILLUMINA" -K 100000000 -M -t 30 $REF \
    ${file} | \
    samtools view -@30 -bS | samtools sort -@30 -m 4G -o $RawData/bam/$dir/${sample/_S*_R1_001.fastq.gz/.bam} - 
    samtools index -@30 $RawData/bam/$dir/${sample/_S*_R1_001.fastq.gz/.bam}
    done

# Remove the reference from memory (only run this once)
bwa shm -d