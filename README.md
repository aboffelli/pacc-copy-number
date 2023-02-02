# PACCs scWGS analysis

## Software installation

All installations are performed in the ~/bin/ directory.

&nbsp;

### Installing GATK

---

Clone from GitHub

```sh
git clone https://github.com/broadinstitute/gatk.git
# The folder contains a stand alone script, so no further installation is necessary.
# Create a link of the script in the bin folder to have the software in PATH. We need to rename the gatk folder first.
mv gatk/ Gatk/; ln -s ~/bin/Gatk/gatk gatk
```

### Instaling Eagle

---

Download tarball

```sh
curl -l https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz -o eagle.tar.gz
tar -xzvf eagle.tar.gz; rm eagle.tar.gz
# Eagle also has a standalone script, so no installation is required, but we can link the script to have it in PATH
ln -s ~/bin/Eagle_v2.4.1/eagle eagle
```

### Installing CHISEL

---

Clone from GitHub and use their installation file

For some reason this works and installing it apart in conda does not...?!

```sh
git clone https://github.com/raphael-group/chisel
cd chisel; bash install_full.sh

# We can test CHISEL with the provided demo
conda activate /home/boffelli/bin/chisel/conda/envs/chisel
cd demos/complete/; bash demo-complete.sh
# A successful run will generate 17 plots and a log file in plots/
```

&nbsp;

## Data preparation

Since we have gotten the BAM files from another lab, we will assume this are raw BAM files.  

### Download reference genome

```sh
# Fasta
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz

# Sort, change nomenclature of chrX and Y index and build a dict from the reference genome
seqkit sort -nN hg38.fa.gz -o hg38_sorted.fa.gz
gunzip hg38_sorted.fa.gz &&
sed -i -e 's/chrX/chr23/' -e 's/chrY/chr24/' hg38_sorted.fa &&
bgzip hg38_sorted.fa &&
samtools faidx hg38_sorted.fa.gz &&
bwa index hg38_sorted_new.fa.gz

bowtie2-build hg38_sorted.fa.gz hg38_sorted
samtools dict hg38_sorted.fa.gz > hg38_sorted.fa.dict

# Download VCF and index
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/20190312_biallelic_SNV_and_INDEL_MANIFEST.txt

cut -f1 20190312_biallelic_SNV_and_INDEL_MANIFEST.txt | sed 's/.\///' | while read f; do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/${f}; done

# Change reference to new nomenclature and convert to bcf
for f in *phased.vcf.gz; do zcat $f | sed -e 's/ID=\([[:digit:]]\)/ID=chr\1/;s/^\([[:digit:]]\)/chr\1/' > ${f/.vcf.gz/_new.vcf}; bcftools view -Ob ${f/.vcf.gz/_new.vcf} > ${f/.vcf/.bcf}; tabix ${f/.vcf.gz/_new.bcf}; rm ${f}*; done

# For already existing new.vcf.gz
for f in *phased.vcf.gz; do bcftools view -Ob ${f} > ${f/.vcf.gz/.bcf}; tabix ${f/.vcf.gz/.bcf}; rm ${f}*; done
```

### List of controls, paccs, and progeny

```sh
# Isolate the control and paccs samples
cat Plate_sample_info.tab | grep "CTR" | cut -f3-5 | while read line; do sample=$(echo ${line} | cut -d ' ' -f1); s=$(echo ${line} | cut -d ' '  -f2); e=$(echo ${line} | cut -d ' ' -f3); for i in $(seq ${s} ${e}); do echo ${sample}/${sample}_${i}; done; done > controls.txt

cat Plate_sample_info.tab | grep "PACC" | cut -f3-5 | while read line; do sample=$(echo ${line} | cut -d ' ' -f1); s=$(echo ${line} | cut -d ' '  -f2); e=$(echo ${line} | cut -d ' ' -f3); for i in $(seq ${s} ${e}); do echo ${sample}/${sample}_${i}; done; done > paccs.txt

cat Plate_sample_info.tab | grep "Progeny" | cut -f3-5 | while read line; do sample=$(echo ${line} | cut -d ' ' -f1); s=$(echo ${line} | cut -d ' '  -f2); e=$(echo ${line} | cut -d ' ' -f3); for i in $(seq ${s} ${e}); do echo ${sample}/${sample}_${i}; done; done > progeny.txt


# fix the numbers
perl -pi -e 's/_(\d)$/_00\1/' controls.txt; perl -pi -e 's/_(\d\d)$/_0\1/' controls.txt

perl -pi -e 's/_(\d)$/_00\1/' paccs.txt; perl -pi -e 's/_(\d\d)$/_0\1/' paccs.txt

perl -pi -e 's/_(\d)$/_00\1/' progeny.txt; perl -pi -e 's/_(\d\d)$/_0\1/' progeny.txt

```

### Align fastq files

```sh
head -2 controls.txt | while read sample; do file=fastq/${sample}*.fastq.gz; bowtie2-align-s --wrapper basic-0 -x RefGenome/hg38_sorted -p 8  --rg-id ${sample#*/} --rg SM:${sample#*/} --rg LB:${sample#*/} --rg PL:Illumina -U ${file} -S newbam/${sample}.sam; done
```

### Convert bam files chromosome notation

```sh
# for f in bam/*/*.bam; do
# new=${f/bam/NewNotationBam};
# samtools view -h $f | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' -e 's/chr\*/\*/' -e 's/chrMT/chrM/' | samtools view -bS - > ${new/%.bam/_new_notation.bam}; 
# done

#  Convert sam to bam and sort.
samtools view -h new3_test.bam | grep -vP "chr\S+?_|chrM" | sed -e 's/chrX/chr23/' -e 's/chrY/chr24/' | samtools view -bo new3_test_new.bam
for f in *.sam; do samtools view -S -b ${f} | samtools sort -o ${f/.sam/_sorted.bam}; done


# Mark  duplicates
for f in *.bam; do java -jar /home/boffelli/bin/picard/build/libs/picard.jar MarkDuplicates -I ${f} -O ${f/.bam/_marked_dup.bam} -M ${f/.bam/_marked_dup_metrics.txt} 2>&1 | tee ${f/.bam/_marked_dup.log}; done

# Index the new bam files
for f in *.bam; do samtools index $f; done
```

---

## Running GATK

```sh
# Germline SNVs ?
head -2 controls.txt | while read sample; do gatk --java-options "-Xmx4g" HaplotypeCaller -R RefGenome/hg38_sorted.fa.gz -I newbam/${sample}*_dup.bam -O GerVCF/${sample}_unfiltered.vcf --native-pair-hmm-threads 8 2>&1 | tee GerVCF/${sample}_unfiltered.log; done


# Run GATK
head -2 controls.txt | while read sample; do gatk Mutect2 -R RefGenome/hg38_sorted.fa.gz -I newbam/${sample}*_dup.bam -O SomVCF/${sample}_unfiltered.vcf 2>&1 | tee SomVCF/${sample}_unfiltered.log; done
head -5 forpaccs.txt | while read sample; do gatk Mutect2 -R RefGenome/hg38.fa.gz -I NewNotationBam/${sample}*.bam -O SomVCF/${sample}.vcf; done

# Filter the VCFs
for f in SomVCF/*/*unfiltered.vcf; do gatk FilterMutectCalls -R RefGenome/hg38_sorted.fa.gz -V ${f} -O ${f/unfiltered/filtered} 2>&1 | tee SomVCF/${sample}_filtered.log; done

```

## Transforming vcf in bcf

```sh
# Convert vcf to bcf and index.
for f in GerVCF/*/*_excesshet.vcf; do bcftools view -f PASS -Ob ${f} > ${f/_excesshet*/.bcf} && tabix ${f/_excesshet*/.bcf}; done

# Split multiallele lines (Not needed if we're filtering when converting to bcf)
for f in SomVCF/*/*.bcf; do bcftools norm -m -any $f > ${f/.bcf/_no_multi.bcf}; done

# Split by chromosome (Maybe not needed? If there is no Y chromosome)
bcftools index -s EH220304_3_096no_multiallele.bcf.bgz | cut -f 1 | while read C; do bcftools view -O z -o split.${C}.vcf.gz EH220304_3_096no_multiallele.bcf.bgz "${C}" ; done
```

## Phasing the controls VCFs

```sh
# Run Eagle for all chromosomes in the control samples. 
f="calls.bcf"
bcftools view -H ${f} | cut -f1 |\
    uniq | grep -vP "chr[24|M]|_" |\
    sed 's/chr//' | while read chr; do
        eagle \
        --vcfRef RefGenome/1000Genomes/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_new.vcf.gz \
        --vcfTarget ${f} \
        --geneticMapFile /home/boffelli/bin/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
        --outPrefix calls_chr${chr}_phased \
        --chrom ${chr} \
        --numThreads 12 2>&1 | tee calls_chr${chr}_phased.log;
    done &&
    bcftools merge calls_*_phased.vcf.gz > f{.bcf/_phased.bcf}.vcf

# Concatenate the files
cat controls.txt | while read sample; do
bcftools concat Phased/${sample}_chr*_phased.vcf.gz | bcftools sort > Phased/${sample}_phased.vcf

```

## Prepare phased tsv

```sh
for f in Phased/*/*_phased.vcf; do grep -oP "^chr\d+\t\d+|\d\|\d(?=:)" $f | perl -0pe 's/\n(?!(chr|$))/\t/g' > ${f/phased.vcf/phases.tsv}; done
```

## Run Chisel

```sh
export TUM="newbam/EH220309_5/EH220309_5_002_sorted_marked_dup_new.bam"
export REF="RefGenome/hg38_sorted_new.fa.gz"
export PHA="Phased/EH220309_5/EH220309_5_002_phases.tsv"
chisel_nonormal -t ${TUM} -r ${REF} -l ${PHA} --seed 12
```
