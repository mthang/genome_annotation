#!/bin/bash

#PBS -P your_poject_code
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N MakeBed

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

export PATH=$PATH:/software/bedops/bin

DATA_DIR=/genome_maize/${GENOME}/EVM/

#GDATA
OUTPUT_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/exonic_regions

mkdir -p ${OUTPUT_DIR}

gff2bed < ${DATA_DIR}/EVM.all.gff3 > ${OUTPUT_DIR}/evm_exon.bed

awk '{print $1"\t"$4"\t"$2"\t"$3}' ${OUTPUT_DIR}/evm_exon.bed | uniq > ${OUTPUT_DIR}/evm_exon_formatted.bed
