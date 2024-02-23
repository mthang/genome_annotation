#!/bin/bash

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N EVM_COMBINE

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
#export SINGULARITY_BIND="/data/folder"

SING_IMAGE_DIR=/singularity

EVM_SIF=${SING_IMAGE_DIR}/evm-1.1.1.sif

cd ${SINGULARITY_BINDPATH}/${GENOME}/EVM

GENOME=${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked

singularity exec ${EVM_SIF} /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

singularity exec ${EVM_SIF} /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome ${GENOME}

find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
