#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -l ncpus=20
#PBS -l mem=100GB
#PBS -l walltime=2:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Kallisto

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

export SINGULARITY_BINDPATH="/path/to/be/mounted/fastq_files"
export SINGULARITY_BIND="/data/folder"

GENOME=Zm-Il14H

KALLISTO=/software/kallisto/kallisto

INDEX=pasa_wo_pseudogenes.idx

${KALLISTO} index -i ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_wo_pseudogenes.idx ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_wo_pseudogenes.fa

${KALLISTO} quant --threads=20 -i ${SINGULARITY_BIND}/${GENOME}/kallisto/${INDEX} -o ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant ${SINGULARITY_BINDPATH}/b73_merged_1.fq.gz ${SINGULARITY_BINDPATH}/b73_merged_2.fq.gz
