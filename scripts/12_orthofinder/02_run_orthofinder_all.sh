#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -l ncpus=30
#PBS -l mem=180GB
#PBS -l walltime=48:00:0
#PBS -l jobfs=1GB
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N OrthoFinder

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

export SINGULARITY_BIND="/path/to/be/mounted"

DATA_DIR=${SINGULARITY_BIND}/data/orthofinder_all

#memory exceeded issue when the number species is growing
singularity exec ${SINGULARITY_BIND}/singularity/OrthoFinder-2.5.4.sif orthofinder -t 30 -f ${DATA_DIR}

