#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=0:10:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Expression

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

GENOME=Zm-Il14H

GFF3_DIR=/genome_maize/

ID_FILE=kallisto/pasa_quant/id_no_expression.txt

awk '$5==0 {print}' ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/abundance.tsv > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/cds_no_expression.tsv

awk '$5!=0 {print}' ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/abundance.tsv > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/cds_w_expression.tsv

cut -f1 ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/cds_no_expression.tsv | sort | uniq | grep -v "target_id" > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/id_no_expression.txt
