#!/bin/bash -l

#SBATCH -A snic2022-22-1233
#SBATCH -p core -n 15 
#SBATCH -J chisel_progeny
#SBATCH -t 06:00:00
#SBATCH --mail-user=arthur.boffelli_castro@med.lu.se
#SBATCH --mail-type=ALL

module load conda
conda activate chisel

export RawData="/home/boffelli/mystorage/Data"
export Running="/home/boffelli/mystorage/running"
export REF="$RawData/RefGenome/hg19.fa"


# Run chisel

chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23" \
	--reference $REF\
	--minreads 30000 --size 10Mb --blocksize 0 --seed 12 \
	--simcov 0.02 --nophasecorr --missingsnps 5,5 --maxploidy 5\
	--tumor $Running/progeny.barcodedcells.q20.bam --listphased $Runnning/calls_phased_hg19.bed 

