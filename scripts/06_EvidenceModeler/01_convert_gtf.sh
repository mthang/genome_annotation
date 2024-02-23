#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=2048MB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Il14HGTF2GFF3

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
BRAKER_DIR=braker_maize

export SINGULARITY_BINDPATH=${BRAKE_DIR}/${GENOME}/annotation_skipOptimize
export PATH=$PATH:/software/gt-1.6.2-Linux_x86_64-64bit-complete/bin

GENEMARK_DIR=/${BRAKER_DIR}/${GENOME}/annotation_skipOptimize/GeneMark-ET

GENEMARK_GTF_GFF=genemark_gtf2gff3.pl

perl ${GENEMARK_GTF_GFF} ${GENEMARK_DIR}/genemark.gtf > ${GENEMARK_DIR}/genemark.gff3

EVM_SIF=/singularity/evm-1.1.1.sif

singularity exec ${EVM_SIF} /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl ${SINGULARITY_BINDPATH}/augustus.hints.gff3 > ${SINGULARITY_BINDPATH}/augustus.hints.evm.gff3
