#!/bin/bash

cl=("HCC" "S786")

# Run the no-normal option of chisel.
for file in "${cl[@]}"; do
    chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23" \
        --reference $REF\
        --minreads 30000 --size 10Mb --seed 12 --simcov 0.02 --maxploidy 5\
        --tumor $Running/$file.barcoded.q15.bam  --listphased $Running/calls_phased_hg19.bed
