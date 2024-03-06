#!/bin/bash

####
#
# This script is used to download the b73 RNAseq FASTQ files
#
####

FTP=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR377

# gff https://maizegdb.org/NAM_project
#https://nam-genomes.github.io/
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8628
#https://www.ebi.ac.uk/ena/browser/view/ERR3773807
# The B73 RNAseq data can be found in the link below
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8628/sdrf

for i in {3773807..3773827}
do
   echo ${i}
   lastdigit=${i: -1}
   echo $lastdigit
   curl -o ERR${i}_1.fastq.gz ${FTP}/00${lastdigit}/ERR${i}/ERR${i}_1.fastq.gz
   curl -o ERR${i}_2.fastq.gz ${FTP}/00${lastdigit}/ERR${i}/ERR${i}_2.fastq.gz
done
