#!/bin/bash

WORKDIR=/path/to/genome/folder


while IFS="" read -r genome || [ -n "$genome" ]
do
    echo $genome
    #mkdir $genome
    cd $genome
   # wget -O ${genome}.genome.tar.gz https://ricerc.sicau.edu.cn/RiceRC/download/downloadFile?name=${genome}.genome.tar.gz
   # wget -O ${genome}.IGDBv1.Allset.gff.tar.gz https://ricerc.sicau.edu.cn/RiceRC/download/downloadFile?name=${genome}.IGDBv1.Allset.gff.tar.gz
   # tar -xzvf ${genome}.genome.tar.gz
   # mv ${genome}.genome ${genome}.genome.fa
    tar -xzvf ${genome}.IGDBv1.Allset.gff.tar.gz
    mv ${genome}/*.gff .
    rm -rf ${genome}/
    cd ..
done < RiceGenome
