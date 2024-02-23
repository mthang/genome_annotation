#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -l ncpus=4
#PBS -l mem=10GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N getFasta

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
export SINGULARITY_BIND="/data/folder"

GENOME=Zm-Il14H
GFF=pasa_pseudogenes_filtered_prekallisto.gff3

mkdir -p ${SINGULARITY_BIND}/${GENOME}/kallisto

# PASA (UTR) after pseudogenes removed
singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/gff3_file_to_proteins.pl ${SINGULARITY_BIND}/${GENOME}/gff/final/${GFF} ${SINGULARITY_BIND}/${GENOME}/${GENOME}.genome.fa CDS > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_wo_pseudogenes.fa
