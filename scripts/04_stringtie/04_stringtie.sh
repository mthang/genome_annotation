#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=5
#PBS -l mem=30GB
#PBS -l walltime=5:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Zm-Il14H_STRT

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

STRINGTIE_GFF_DIR=/original/gff/file/of/Zm-Il14H
# path to stringtie software
STRINGTIE=/bin/stringtie

HISAT_BAM=${STRINGTIE_GFF_DIR}/${GENOME}/hisat/b73_sorted_all.bam

#Temp directory
export TMPDIR=/scratch/tmp

#https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/

mkdir -p ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie

${STRINGTIE} ${HISAT_BAM} -l merged -p 5 -G ${STRINGTIE_GFF_DIR}/${GENOME}/${GENOME}.gff3 -o ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/merged_stringtie.gtf

${STRINGTIE} --merge -p 5 -G ${STRINGTIE_GFF_DIR}/${GENOME}/${GENOME}.gff3 -o ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/final_merged_stringtie.gtf ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/merged_stringtie.gtf

