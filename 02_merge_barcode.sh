#!/bin/bash

# Remember to activate chisel env
conda activate chisel

cl=("HCC" "S786")

# merge bam files and barcode.
for dir in "${cl[@]}"; do
    chisel_prep $RawData/bam/$dir/*.bam \
    -o $Running/$dir.barcoded.bam 

    # Filter read quality.
    samtools view -@ 30 -q 20 $Running/$dir.barcoded.q20.bam $Running/$dir.barcoded.bam
done