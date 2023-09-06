#!/bin/bash

# First, create a map file to sort the plots. 
map="HCC"
info_file="/home/boffelli/paccs/running/AlignmentQC/${map}.barcoded.info.tsv"

# Use the Plate sample info file to get the cells that belong to each cell line.
cat $HOME/raw/Plate_sample_info.tab | grep ${map} | cut -f3-6 | \
while read line; do 
    sample=$(echo ${line} | cut -d ' ' -f1); 
    s=$(echo ${line} | cut -d ' '  -f2); 
    e=$(echo ${line} | cut -d ' ' -f3); 
    rep=$(echo ${line} | cut -d ' ' -f4); 
    for i in $(seq ${s} ${e}); do 
        echo -e ${sample}_${i}"\t"$rep; 
        done; 
done | sort -k 2 | perl -p -e 's/_(\d\t)/_00\1/' | perl -p -e 's/_(\d\d\t)/_0\1/' | grep -v "Progeny" |\

while read line; do 
    sample=$(echo $line | grep -oP "EH\d+_\d_\d+"); 
    bar=$(echo $line | grep "$sample" ${map}.barcoded.info.tsv | cut -f2); 
    clone=$(echo $line | cut -d ' ' -f2); 
    echo -e "$bar\t$clone" | \
    # Uncomment the next line according to the cell line you are using
    sed -e 's/HCC_CTR_R1/1\tHCC_CTR_R1/' -e 's/HCC_CTR_R2/2\tS786_O_CTR_R2/' -e 's/HCC_CTR_R3/3\tS786_O_CTR_R3/' -e 's/HCC_PACCs_R1/4\tHCC_PACCs_R1/' -e 's/HCC_PACCs_R2/5\tHCC_PACCs_R2/' -e 's/HCC_PACCs_R3/6\tHCC_PACCs_R3/';
    #sed -e 's/S786_O_CTR_R1/1\tS786_O_CTR_R1/' -e 's/S786_O_CTR_R2/2\tS786_O_CTR_R2/' -e 's/S786_O_CTR_R3/3\tS786_O_CTR_R3/' -e 's/S786_O_PACCs_R1/4\tS786_O_PACCs_R1/' -e 's/S786_O_PACCs_R2/5\tS786_O_PACCs_R2/' -e 's/S786_O_PACCs_R3/6\tS786_O_PACCs_R3/'; 
done | grep -vP "^\t" > ${map}.map

# Remove problematic parts from the calls file
# Inside the calls directory created by chisel
cd calls
bedtools intersect -header -a calls.tsv -b ../keep-hg19_c23.bed -wa | uniq > calls.tsv.new
mv calls.tsv.new calls.tsv

# Back to running directory, rerun the chisel plot.
cd ..
chisel_plotting -m ${map}.map -figformat pdf