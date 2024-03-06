#!/bin/bash


export DATA_DIR="/path/to/fastq/folder"

cd ${DATA_DIR}


while IFS="" read -r line || [ -n "$line" ]
do
    #echo $line
    #mkdir -p repeat/$genome
    #SPECIES=`echo $line | awk '{print $1}' | cut -d"_" -f1`
    SPECIES=`echo $line | awk '{print $1}' | sed 's/_data//g'`
    R1=`echo $line | awk '{print $2}'`
    R2=`echo $line | awk '{print $3}'`
    R1_name=`basename ${R1}`
    R2_name=`basename ${R2}`
    #echo -e "${SPECIES} ${R1} ${R2}"
    echo -e "${SPECIES}_${R1_name}"
    wget -O ${SPECIES}_${R1_name} ${R1}
    echo -e "${SPECIES}_${R2_name}"
    wget -O ${SPECIES}_${R2_name} ${R2}
done < RNASEQ_DATA.txt
