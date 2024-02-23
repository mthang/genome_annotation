#!/bin/bash

#PBS -P kw68
#PBS -l ncpus=1
#PBS -l mem=5GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N 7BmakeFolder

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

DATA_DIR=/genome_maize/${GENOME}/pseudogenes/EVM

CHR_DIR=${DATA_DIR}/chr
EXON_DIR=${DATA_DIR}/exonic_regions

PROTEIN=/genome_maize/${GENOME}/EVM/evm_prot.fa

cd ${DATA_DIR}

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr}
    ls ${CHR_DIR}/${chr}.fa
    ls ${CHR_DIR}/${chr}_rm.fa
    ls ${EXON_DIR}/${chr}.bed
    ls ${PROTEIN}

    mkdir -p ${chr}/input/pep
    mkdir -p ${chr}/input/dna
    mkdir -p ${chr}/input/mysql

        cp ${CHR_DIR}/${chr}.fa ${chr}/input/dna/
    cp ${CHR_DIR}/${chr}_rm.fa ${chr}/input/dna/
    cp ${PROTEIN} ${chr}/input/pep/pep.fa
    cp ${EXON_DIR}/${chr}.bed ${chr}/input/mysql/${chr}_exLocs
done < ${DATA_DIR}/chr/chr_id.txt

    
