#!/bin/bash -l

#SBATCH -A snic2022-22-1233
#SBATCH -p core
#SBATCH -n 15
#SBATCH -t 06:00:00
#SBATCH -J alignment
#SBATCH --mail-user=arthur.boffelli_castro@med.lu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bwa
module load samtools
#source /home/boffelli/bin/chisel/conda/bin/activate chisel

# BWA mem alignment.
export RawData="/home/boffelli/mystorage/private/Data"
export Running="/home/boffelli/mystorage/private/running"
#for file in $RawData/controls/*.fastq.gz; do 
#    sample=${file##*/}
#    bwa mem -t 30 $RawData/RefGenome/hg38_sorted_new.fa.gz \
#    ${file} | \
#    samtools sort -@30 -o $Running/bam/controls/${sample/_S*_R1_001.fastq.gz/.bam}
#    done
#for file in $RawData/progeny/*.fastq.gz; do 
#    sample=${file##*/}
#    bwa mem -t 30 $RawData/RefGenome/hg38_sorted_new.fa.gz \
#    ${file} | \
#    samtools sort -@30 -o $Running/bam/progeny/${sample/_S*_R1_001.fastq.gz/.bam}
#    done
for file in $RawData/paccs/*.fastq.gz; do 
    sample=${file##*/}
    bwa mem -t 30 $RawData/RefGenome/hg38_sorted_new.fa.gz \
    ${file} | \
    samtools sort -@30 -o $Running/bam/paccs/${sample/_S*_R1_001.fastq.gz/.bam}
    done


## merge bam files and barcode.
#chisel_prep ../Data/bam/controls/*.bam \
#	-o ../Data/merged.barcoded.controls.bam 
#
# Variant call only for controls. 
#bcftools mpileup --ignore-RG -Ob \
#	-f $HOME/pacc_raw/RefGenome/hg38_sorted_new.fa.gz merged.controls.barcoded.bam |\
# #	bcftools call -mv --threads 30 | bcftools filter -i 'QUAL>30' -Ob -o calls.bcf
# 
#f="calls.bcf"
# 
# #bcftools index ${f}
# 
# Phase the variants.
#bcftools view -H ${f} | cut -f1 |\
#    uniq | grep -vP "chr24|chrM|_"|\
#    sed 's/chr//' | while read chr; do
#        eagle \
#        --vcfRef $HOME/pacc_raw/RefGenome/1000Genomes/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_new.bcf \
#        --vcfTarget ${f} \
#        --geneticMapFile /home/boffelli/bin/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
#        --outPrefix calls_chr${chr}_phased \
#        --chrom ${chr} \
#        --numThreads 30 2>&1 | tee calls_chr${chr}_phased.log
#        done
#bcftools concat calls_*_phased.vcf.gz | bcftools sort > ${f/.bcf/_phased.vcf} && rm calls_*_phased.vcf.gz
#grep -oP "^chr\d+\t\d+|\d\|\d(?=:)" calls_phased.vcf | perl -0pe 's/\n(?!(chr|$))/\t/g' > phases.tsv
# Run chisel
#export TUM="merged.progeny.barcoded.bam"
#export REF="$HOME/pacc_raw/RefGenome/hg38_sorted_new.fa.gz"
#export PHA="phases.tsv" 
#
#chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22" \
#	--reference $REF\
#	--minreads 30000 --size 10Mb --blocksize 0 --seed 12 \
#    --simcov 0.02 --nophasecorr --missingsnps 5,5 --maxploidy 5\
#	--tumor $TUM --listphased $PHA 

