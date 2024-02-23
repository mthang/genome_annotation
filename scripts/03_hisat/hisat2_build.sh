#!/bin/bash


##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=20
#PBS -l mem=15GB
#PBS -l walltime=00:20:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Zm-Ky21

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

export SINGULARITY_BINDPATH="/path/to/be/mounted"

GENOME=Zm-Il14H

singularity exec ${SINGULARITY_BINDPATH}/hisat2.sif hisat2-build -p 20 ${SINGULARITY_BINDPATH}/${GENOME}/${GENOME}.genome.fa ${SINGULARITY_BINDPATH}/${GENOME}/${GENOME}.genome.fa
