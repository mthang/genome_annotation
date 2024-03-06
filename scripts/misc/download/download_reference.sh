#!/bin/bash

URL=https://download.maizegdb.org
#https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz

while IFS="" read -r line || [ -n "$line" ]
do
    echo $line
    folder=`basename $line .fa.gz`
    species=`echo $folder | sed 's/-REFERENCE-NAM-[15].0//g'`
    echo -e "$folder $species"
    mkdir $species
    curl -o ${species}/$species.fa.gz "${URL}/${folder}/${line}"
    #wget -O ${species}/${species}.gff3.gz --recursive --level=1 --no-parent --no-directories --accept ${folder}_*.1.gff3.gz ${URL}/${folder}/
    wget --recursive --level=1 --no-parent --no-directories --accept ${folder}_*.1.gff3.gz ${URL}/${folder}/
    mv *.gff3.gz ${species}/${species}.gff3.gz
    gunzip ${species}/*.gz
done < maize.txt
