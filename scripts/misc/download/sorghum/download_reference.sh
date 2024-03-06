#!/bin/bash


#####
#
# This script is used to download sorghum genome fasta files 
#
#####

URL=https://ftp.cngb.org/pub/CNSA/data3/CNP0001440

## 13 genomes
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314252/CNA0019255/353.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314253/CNA0019257/IS12661.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314254/CNA0019259/IS3614-3.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314255/CNA0019261/IS929.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314256/CNA0019256/AusTRCF317961.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314257/CNA0019258/IS19953.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314258/CNA0019260/IS8525.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314259/CNA0019262/Ji2731.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314260/CNA0019263/PI525695.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314261/CNA0019264/PI532566.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314262/CNA0019265/PI536008.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314263/CNA0019266/R931945-2-2.genome.fa.gz
# https://ftp.cngb.org/pub/CNSA/data3/CNP0001440/CNS0314264/CNA0019267/S369-1.genome.fa.gz

## other three (Rio, BTx) genomes are availble in the link below
# 1) go to https://www.sorghumbase.org/post/sorghum-bicolor-btx623
# 2) click on "GENOME"
# 3) select genome of interest (i.e BTx623, Rio and TX430)
# 4) use this url https://ensembl.sorghumbase.org/Sorghum_rio/Info/Index?db=core to download both fasta and gff3 file
# 5) repeat the steps above for different genomes
# Rio - Sorghum_rio.JGI-v2.0.dna.toplevel.fa
# BTx623 - Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
# TX430 - Sorghum_tx430nano.Corteva_Sorghum_ONT_TX430_1.0.dna.toplevel.fa

FILES=(CNS0314252/CNA0019255/353.genome.fa.gz 
CNS0314253/CNA0019257/IS12661.genome.fa.gz 
CNS0314254/CNA0019259/IS3614-3.genome.fa.gz
CNS0314255/CNA0019261/IS929.genome.fa.gz
CNS0314256/CNA0019256/AusTRCF317961.genome.fa.gz
CNS0314257/CNA0019258/IS19953.genome.fa.gz
CNS0314258/CNA0019260/IS8525.genome.fa.gz
CNS0314259/CNA0019262/Ji2731.genome.fa.gz
CNS0314260/CNA0019263/PI525695.genome.fa.gz
CNS0314261/CNA0019264/PI532566.genome.fa.gz
CNS0314262/CNA0019265/PI536008.genome.fa.gz
CNS0314263/CNA0019266/R931945-2-2.genome.fa.gz
CNS0314264/CNA0019267/S369-1.genome.fa.gz)

# reverse_end=`echo ${FASTQ_DIR}/${p} | sed 's/_1/_2/g'`

for f in ${FILES[@]}
do
   if [[ "${f}" != *.fa.gz && "${f}" != *.gfa && "${f}" != *.txt ]];then
       echo ${URL}/${f}
       wget -nH --cut-dirs=3 -x ${URL}/${f}
       if [[ "${f}" == *_1.fq.gz ]]; then
          reverse_end=`echo ${f} | sed 's/_1/_2/g'`
          echo ${URL}/${reverse_end}
          wget -nH --cut-dirs=3 -x ${URL}/${reverse_end}
       fi
   else
       echo ${URL}/${f}
       echo "downloading misc"
       wget -nH --cut-dirs=3 -x ${URL}/${f}
   fi
done
