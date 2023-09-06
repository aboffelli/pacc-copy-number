# PACCs scWGS analysis

## Software installation

All installations are performed in the ~/bin/ directory.



### Instaling Eagle

```sh
curl -l https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz -o eagle.tar.gz
tar -xzvf eagle.tar.gz; rm eagle.tar.gz
# Eagle also has a standalone script, so no installation is required, but we can link the script to have it in PATH
ln -s ~/bin/Eagle_v2.4.1/eagle eagle
```

### Installing CHISEL

```sh
git clone https://github.com/raphael-group/chisel
cd chisel; bash install_full.sh

# We can test CHISEL with the provided demo
conda activate /home/boffelli/bin/chisel/conda/envs/chisel
cd demos/complete/; bash demo-complete.sh
# A successful run will generate 17 plots and a log file in plots/
```

&nbsp;

### Download reference genome

All the reference genome and vcf files used can be found in this ftp link:  
**ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37**

The only changes made to the reference genome was changing the nomenclature of chrX to chr23 and chrY to chr24, for better usage in the software. After changing the nomenclature with, e.g. ```sed -i -e 's/chrX/chr23/' -e 's/chrY/chr24/' hg19.fa > hg19_new.fa``` you need to index and build the dictionary again.
Using for example, ```samtools faidx hg19_new.fa && samtools dict hg19_new.fa > hg19_new.dict``` 

&nbsp;

## 01. Align the fastq files

The alignment of FASTQ files to the reference genome is in ***01_read_alignment.sh*** script. This script will align all samples from control, PACC and progeny.

&nbsp;

## 02. Merge the bam files and barcode the cells

Use Chisel-prep command to merge all cells from each cell line together and barcode each sample. Use script ***02_merge_barcode.sh***.

&nbsp;

## 03. Phasing the controls VCFs. 

We only need the haplotype phasing in the SNPs of the control cells. The script ***03_vc_phasing.sh*** contains the code for this section.

&nbsp;

## 04. Chisel no-normal

Run the Chisel no-normal option, since we don't have normal cells to compare with. Described in the script ***04_chisel_nonormal.sh***.

&nbsp;

## 05.Plotting

Chisel will already generate plots after running the previous step, however we can sort the cells and remove unwanted regions from the plot. For that use the script ***05_chisel_plot.sh***

The plots will be located in the ***plots*** folder, created by Chisel. Our plot of interest is the ***rdrs.pdf***.
