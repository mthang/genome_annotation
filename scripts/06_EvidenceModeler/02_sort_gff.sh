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
#PBS -N SortGff3

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

PROJECT_DIR=/genome_maize/${GENOME}

export PATH=$PATH:/software/gt-1.6.2-Linux_x86_64-64bit-complete/bin

BRAKER_DIR=/braker_maize/${GENOME}/annotation_skipOptimize

GENEMARK_GFF=${BRAKER_DIR}/GeneMark-ET/genemark.gff3

AUGUSTUS_GFF=${BRAKER_DIR}/augustus.hints.evm.gff3

SNAP_GFF=${PROJECT_DIR}/snap/snap.gff

FGENESH_GFF=${PROJECT_DIR}/fgenesh/${GENOME}.gff3

gt gff3 -sort -tidy -retainids ${GENEMARK_GFF} > ${BRAKER_DIR}/GeneMark-ET/genemark_sorted.gff3

gt gff3 -sort -tidy -retainids ${AUGUSTUS_GFF} > ${BRAKER_DIR}/augustus.hints.evm.sorted.gff3

gt gff3 -sort -tidy -retainids ${FGENESH_GFF} > ${PROJECT_DIR}/fgenesh/${GENOME}_sorted.gff3

gt gff3 -sort -tidy -retainids ${SNAP_GFF} > ${PROJECT_DIR}/snap/snap_sorted.gff3
