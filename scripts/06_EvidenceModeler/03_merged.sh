#!/bin/bash


##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=3GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N MergeGFF3

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

#######
# DO NOT NEED TO RUN THIS SCRIPT AGAIN IF THE MERGED.GFF3 ALREADY CREATED #####
#######
GENOME=Zm-Il14H

PROJECT_DIR1=/g/data/kw68/genome_maize/${GENOME}
PROJECT_DIR2=/scratch/kw68/wt5249/

OUTPUT_DIR=${PROJECT_DIR1}/gff/merged

mkdir -p ${OUTPUT_DIR}

AUGUSTUS_GFF=${PROJECT_DIR2}/braker_maize/${GENOME}/annotation_skipOptimize/augustus.hints.evm.gff3
GENEMARK_GFF=${PROJECT_DIR2}/braker_maize/${GENOME}/annotation_skipOptimize/GeneMark-ET/genemark_sorted.gff3

#PASA_GFF=${PROJECT_DIR2}/pasa/${GENOME}/sample_data/${GENOME}.sqlite.pasa_assemblies.gff3

SNAP_GFF=${PROJECT_DIR1}/snap/snap_sorted.gff3

#FGENESH_GFF=${PROJECT_DIR1}/fgenesh/${GENOME}_sorted.gff3

grep -v "^# \|^###" ${SNAP_GFF} > ${OUTPUT_DIR}/snap_processed.gff

grep -v "^# \|^###" ${GENEMARK_GFF} > ${OUTPUT_DIR}/genemark_processed.gff

# comment out if fgenesh_processed.gff is already processed (i.e sorted)
# grep -v "^# \|^###" ${FGENESH_GFF} > ${OUTPUT_DIR}/${GENOME}_fgenesh_processed.gff

ln -s ${PROJECT_DIR1}/fgenesh/fgenesh_processed.gff ${OUTPUT_DIR}/${GENOME}_fgenesh_processed.gff

cat ${OUTPUT_DIR}/${GENOME}_fgenesh_processed.gff ${AUGUSTUS_GFF} ${OUTPUT_DIR}/genemark_processed.gff ${OUTPUT_DIR}/snap_processed.gff > ${OUTPUT_DIR}/merged_tmp.gff

grep -v "^#" ${OUTPUT_DIR}/merged_tmp.gff > ${OUTPUT_DIR}/merged.gff

rm ${OUTPUT_DIR}/merged_tmp.gff

