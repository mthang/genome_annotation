#!/bin/bash

#PBS -P your_project_code
#PBS -l ncpus=2
#PBS -l mem=3GB
#PBS -l walltime=3:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N POST_EVM_PROTEIN

module load singularity

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

export SINGULARITY_BINDPATH="/path/to/be/mounted"
export SINGULARITY_BIND="/data/folder/"

SING_IMAGE_DIR=/singularity

GFF3=${SINGULARITY_BINDPATH}/${GENOME}/EVM/EVM.all.gff3

FASTA=${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked

singularity exec ${SING_IMAGE_DIR}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/gff3_file_to_proteins.pl  ${GFF3} ${FASTA} > ${SINGULARITY_BINDPATH}/${GENOME}/EVM/evm_prot.fa
