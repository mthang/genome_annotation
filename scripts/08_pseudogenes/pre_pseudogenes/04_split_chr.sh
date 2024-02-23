#!/bin/bash

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=5GB
#PBS -l walltime=5:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N splitfolder

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

GENOME=Zm-Il14H

DATA_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/exonic_regions
# GENOME
SPECIES=/genome_maize/${GENOME}/${GENOME}.genome.fa
# Genome from repeat masking
SPECIES_RM=/genome_maize/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked
OUTPUT_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/chr

mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

cut -f1 ${DATA_DIR}/evm_exon.bed | sort | uniq > ${OUTPUT_DIR}/chr_id.txt

FASOMERECORDS=/software/faSomeRecords

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr} > tmp.txt
    ${FASOMERECORDS} ${SPECIES} tmp.txt ${chr}.fa
    ${FASOMERECORDS} ${SPECIES_RM} tmp.txt ${chr}_rm.fa
    rm tmp.txt
done < chr_id.txt
