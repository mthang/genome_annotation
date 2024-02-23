#!/bin/bash

#PBS -P your_project_code
#PBS -l ncpus=20
#PBS -l mem=20GB
#PBS -l walltime=25:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N CML333_RUN_EVM

module load singularity
module load parallel/20191022

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

parallel -j 20 --workdir ${SINGULARITY_BINDPATH}/${GENOME}/EVM < commands.list
