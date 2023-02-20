#!/bin/bash -l

#SBATCH -A snic2022-22-1233
#SBATCH -p core
#SBATCH -n 15
#SBATCH -t 7-00:00:00
#SBATCH -J whole_thing
#SBATCH --mail-user=arthur.boffelli_castro@med.lu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bwa
module load samtools
module load bcftools
module load conda

# BWA mem alignment.
export RawData="/home/boffelli/mystorage/Data"
export Running="/home/boffelli/mystorage/running"
export REF="$RawData/RefGenome/hg19.fa"

bwa shm $REF

for file in $RawData/fastq/controls/*.fastq.gz; do 
    sample=${file##*/}
    bwa mem -R "@RG\tID:$sample\tSM:control\tLB:control\tPL:ILLUMINA\tPU:ILLUMINA" -K 100000000 -M -t 30 $REF \
    ${file} | \
    samtools view -@30 -bS | samtools sort -@30 -m 4G -o $RawData/bam/controls/${sample/_S*_R1_001.fastq.gz/.bam} -
    samtools index -@30 $RawData/bam/controls/${sample/_S*_R1_001.fastq.gz/.bam}
    done
for file in $RawData/fastq/paccs/*.fastq.gz; do 
    sample=${file##*/}
    bwa mem -R "@RG\tID:$sample\tSM:pacc\tLB:pacc\tPL:ILLUMINA\tPU:ILLUMINA" -K 100000000 -M -t 30 $REF \
    ${file} | \
    samtools view -@30 -bS | samtools sort -@30 -m 4G -o $RawData/bam/paccs/${sample/_S*_R1_001.fastq.gz/srt.bam} -
    samtools index -@30 $RawData/bam/paccs/${sample/_S*_R1_001.fastq.gz/srt.bam}
    done
#for file in $RawData/fastq/progeny/*.fastq.gz; do 
#    sample=${file##*/}
#    bwa mem -R "@RG\tID:$sample\tSM:progeny\tLB:progeny\tPL:ILLUMINA\tPU:ILLUMINA" -K 100000000 -M -t 30 $REF \
#    ${file} | \
#    samtools view -@30 -bS | samtools sort -@30 -m 4G -o $RawData/bam/progeny/${sample/_S*_R1_001.fastq.gz/.bam} -
#    samtools index -@30 $RawData/bam/progeny/${sample/_S*_R1_001.fastq.gz/.bam}
#    done
bwa shm -d

conda activate chisel

# merge bam files and barcode.
chisel_prep $RawData/bam/controls/*.bam \
	-o $Running/controls.barcoded.bam 

chisel_prep $RawData/bam/paccs/*.bam \
    -o $Running/paccs.barcoded.bam

#chisel_prep $RawData/bam/progeny/*.bam \
#	-o $Running/progeny.barcoded.bam 

# Filter read quality.
samtools view -@ 30 -q 20 $Running/controls.barcoded.q20.bam $Running/controls.barcoded.bam
samtools view -@ 30 -q 20 $Running/paccs.barcoded.q20.bam $Running/paccs.barcoded.bam
#samtools view -@ 30 -q 20 $Running/progeny.barcoded.q20.bam $Running/progeny.barcoded.bam

## Variant call only for controls. 
#bcftools mpileup --ignore-RG -Ob \
#	-f $REF merged.barcoded.controls.bam |\
#	bcftools call -mv --threads 30 | bcftools filter -i 'QUAL>30' -Ob -o calls.bcf
# 
#f="calls.bcf"
# 
#bcftools index ${f}
# 
## Phase the variants.
#bcftools view -H ${f} | cut -f1 |\
#    uniq | grep -vP "chr24|chrM|_"|\
#    sed 's/chr//' | while read chr; do
#        eagle \
#        --vcfRef $RawData/RefGenome/1000Genomes/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_new.bcf \
#        --vcfTarget ${f} \
#        --geneticMapFile $HOME/bin/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
#        --outPrefix calls_chr${chr}_phased \
#        --chrom ${chr} \
#        --numThreads 30 2>&1 | tee calls_chr${chr}_phased.log
#        done
#bcftools concat calls_*_phased.vcf.gz | bcftools sort > ${f/.bcf/_phased.vcf}
#
#grep -oP "^chr\d+\t\d+|\d\|\d(?=:)" calls_phased.vcf | perl -0pe 's/\n(?!(chr|$))/\t/g' | grep -v "1|1" > phases.tsv

# Run chisel

chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23" \
        --reference $REF\
        --minreads 30000 --size 10Mb --seed 12 --simcov 0.02 --maxploidy 5\
        --tumor $Running/controls.barcoded.q15.bam  --listphased $Running/calls_phased_hg19.bed
chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23" \
        --reference $REF\
        --minreads 30000 --size 10Mb --seed 12 --simcov 0.02 --maxploidy 8\
        --tumor $Running/paccs.barcoded.q15.bam  --listphased $Running/calls_phased_hg19.bed
#chisel_nonormal --chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23" \
#        --reference $REF\
#        --minreads 30000 --size 10Mb --seed 12 --simcov 0.02 --maxploidy 5\
#        --tumor $Running/progeny.barcoded.q15.bam  --listphased $Running/calls_phased_hg19.bed

