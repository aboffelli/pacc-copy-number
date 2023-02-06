#!/bin/bash -l

#SBATCH -A snic2022-22-1233
#SBATCH -p core
#SBATCH -n 15
#SBATCH -t 7-00:00:00
#SBATCH -J alignment_variant_phasing
#SBATCH --mail-user=arthur.boffelli_castro@med.lu.se
#SBATCH --chdir=/home/boffelli/mystorage/private
#SBATCH --output=/home/boffelli/mystorage/private/slurm_out
#SBATCH --error=/home/boffelli/mystorage/private/slurm_err

source /home/boffelli/bin/chisel/conda/bin/activate chisel

# New version

# TODO Add bwa mem alignment.
# Test
#head -5 $HOME/pacc_raw/controls.txt | while read sample; do 
#    bwa mem -t 12 $HOME/pacc_raw/RefGenome/hg38_sorted_new.fa.gz \
#    $HOME/pacc_raw/fastq/$sample*.fastq.gz | \
#    samtools sort -@12 -o $HOME/paccs/bam/controls/${sample#*/}.bam 2>&1 | tee $HOME/paccs/bam/controls/${sample#*/}.log
#    done
#head -5 $HOME/pacc_raw/progeny.txt | while read sample; do 
#    bwa mem -t 12 $HOME/pacc_raw/RefGenome/hg38_sorted_new.fa.gz \
#    $HOME/pacc_raw/fastq/$sample*.fastq.gz | \
#    samtools sort -@12 -o $HOME/paccs/bam/progeny/${sample#*/}.bam 2>&1 | tee $HOME/paccs/bam/progeny/${sample#*/}.log
#    done

## merge bam files and barcode.
#chisel_prep ../Data/bam/controls/*.bam \
#	-o ../Data/merged.barcoded.controls.bam 
#
# Variant call only for controls. 
bcftools mpileup --ignore-RG -Ob \
	-f $HOME/pacc_raw/RefGenome/hg38_sorted_new.fa.gz merged.controls.barcoded.bam |\
	bcftools call -mv --threads 30 | bcftools filter -i 'QUAL>30' -Ob -o ../calls.bcf

f="calls.bcf"

bcftools index ${f}

# Phase the variants.
#bcftools view -H ${f} | cut -f1 |\
#    uniq | grep -vP "chr24|chrM|_"|\
#    sed 's/chr//' | while read chr; do
#        eagle \
#        --vcfRef RefGenome/1000Genomes/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_new.bcf \
#        --vcfTarget ${f} \
#        --geneticMapFile /home/boffelli/bin/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
#        --outPrefix calls_chr${chr}_phased \
#        --chrom ${chr} \
#        --numThreads 30 2>&1 | tee calls_chr${chr}_phased.log;
#    done && \
#    bcftools concat calls_*_phased.vcf.gz | bcftools sort > $f{.bcf/_phased.bcf}.vcf
#
## Run chisel
#export TUM="progeny.bam"
#export NOR="controls.bam"
#export REF="RefGenome/hg38_sorted_new.fa.gz"
#export PHA="calls_phased.vcf" 
#
#chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23" \
#	--reference $REF\
#	--minreads 30000 --size 10Mb --blocksize 0 --seed 12 --simcov 0.02 --nophasecorr --missingsnps 5,5\
#	--tumor $TUM --listphased $PHA 
#
