#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N MikadoStats

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
export SINGULARITY_BIND="/path/to/be/mounted"

#Zm-Il14H        967772
GENOME=Zm-Il14H
GENOME_DIR=/genome_maize
PASA_NUM=967772

GFF3=${GENOME}.sqlite.gene_structures_post_PASA_updates.${PASA_NUM}.gff3

singularity exec /singularity/mikado-2.3.3.sif mikado util grep ${SINGULARITY_BIND}/${GENOME}/gff/final/keep.list ${SINGULARITY_BINDPATH}/pasa/${GFF3} ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.gff3

singularity exec /singularity/mikado-2.3.3.sif mikado util stats ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.gff3 ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.stats
