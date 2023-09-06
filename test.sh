#!/bin/bash -l

export REF="/home/boffelli/raw/RefGenome/hg19.fa"

declare -i chr1=4
declare -i chr2=4
declare -i chr19=4
declare -i chr23=2

# Run chisel
for dir in {chr1,chr2,chr19,chr23}; do
cd $dir
chisel_nonormal --chromosomes "chr1" \
        --reference $REF\
        --minreads 30000 --size 10Mb --seed 12 --simcov 0.02 --maxploidy 4\
        --tumor ../parent.barcoded.q20.bam  --listphased ../calls_phased_hg19.bed 
cd ..
done
