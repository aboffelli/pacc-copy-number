#!/bin/bash -l

#SBATCH -A snic2022-22-1233
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 7-00:00:00
#SBATCH -J alignment_variant_phasing
#SBATCH --mail-user=arthur.boffelli_castro@med.lu.se
#SBATCH --chdir=/home/boffelli/mystorage/private
#SBATCH --output=/home/boffelli/mystorage/private/slurm_out
#SBATCH --error=/home/boffelli/mystorage/private/slurm_err

source /home/boffelli/bin/chisel/conda/bin/activate chisel
rm -r _* *.tmp
# New version
# chisel_prep controls/*.fastq.gz \
# 	-r ../RefGenome/hg38_sorted_new.fa.gz \
# 	-o controls.bam \
# 	--rexpname "(.*)_S.*_R[1|2]_001.fastq.*" \
# 	--rexpread ".*_S.*_(R[1|2])_001.fastq.*"
# 
# bcftools mpileup --ignore-RG -Ob \
# 	-f RefGenome/hg38_sorted_new.fa.gz controls.bam |\
#        	bcftools call -mv | bcftools filter -i'QUAL>30' -Ob -o calls.bcf

f="calls.bcf"
# bcftools index ${f}
bcftools view -H ${f} | cut -f1 |\
    uniq | grep -vP "chr24|chrM|_"|\
    sed 's/chr//' | while read chr; do
        eagle \
        --vcfRef RefGenome/1000Genomes/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_new.bcf \
        --vcfTarget ${f} \
        --geneticMapFile /home/boffelli/bin/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
        --outPrefix calls_chr${chr}_phased \
        --chrom ${chr} \
        --numThreads 12 2>&1 | tee calls_chr${chr}_phased.log;
    done && \
    bcftools concat calls_*_phased.vcf.gz | bcftools sort > $f{.bcf/_phased.bcf}.vcf

export TUM="progeny.bam"
export NOR="controls.bam"
export REF="RefGenome/hg38_sorted_new.fa.gz"
export PHA="calls_phased.vcf" 

chisel -t $TUM -n $NOR -r $REF -l $PHA -c "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24"

chisel_nonormal -t $TUM -r $REF -l $PHA --seed 12 -c "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24"
