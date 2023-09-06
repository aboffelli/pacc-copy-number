#!/bin/bash

# Change paths according to your data
GATK_PATH=$HOME/bin/GenomeAnalysisTK.jar
GATK4_PATH=$HOME/gatk-4.0.8.1/gatk
ROOTDIR=$HOME/human_db/GATK_b37/
REF=$ROOTDIR/human_g1k_v37.fasta
DBSNP=$ROOTDIR/dbsnp_138.b37.vcf.gz

# Bam file sample
Tsample=bam_file_name

# Variant calling
$GATK4_PATH HaplotypeCaller -R $REF --dbsnp $DBSNP  -I $Tbam --output out.vcf
bgzip out.vcf
tabix out.vcf.gz

# Select only SNPs and filter quality.
# Since we are only interested in the heterozigous SNPs, we can select only GT="het"
bcftools filter -g 10 out.vcf.gz | bcftools view -v snps  \
    -i 'INFO/DP>=20.0 & MQ>35.0' \
    -i 'GT="het"' - > out.het.vcf

# After that, extract each chromosome data from out.het.vcf file


# Eagle program for haplotype phase
# https://alkesgroup.broadinstitute.org/Eagle/#x1-310005.3.1

# Haplotype phase code
# Run this for every chromosome vcf file.
PANEL_PATH=$HOME/phasing/imputation/1000K_panel

$HOME/phasing/Eagle_v2.4.1/eagle \
--Kpbwt 10000 \
--pbwtIters 10 \
--expectIBDcM 0.2 \
--numThreads 16 \
--vcfRef=$PANEL_PATH/ALL.chr9.phase3_integrated.20130502.genotypes.bcf \
--vcfTarget=chr9.vcf.gz \
--geneticMapFile=$HOME/phasing/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--outPrefix=hete.phased.02cM 2>&1 | tee ref.log

# Combine all chromosomes’ data