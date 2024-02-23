#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=24
#PBS -l mem=12GB
#PBS -l walltime=48:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Zm-Il14H

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

export SINGULARITY_BINDPATH="/path/to/be/mounted/"
export DATA_DIR="/species/folder"

# Species 
SPECIES=Zm-Il14H

cd ${DATA_DIR}/${SPECIES}

# Step 1
singularity exec ${SINGULARITY_BINDPATH}/tetools_repeat.sif BuildDatabase -name ${SPECIES} -engine ncbi ${SPECIES}.fa

#RepeatModeler -database IS8525 -engine ncbi -pa 10
# Step 2
singularity exec ${SINGULARITY_BINDPATH}/tetools_repeat.sif RepeatModeler -engine ncbi -threads 24 -database ${SPECIES}
