#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=2
#PBS -l mem=8GB
#PBS -l walltime=48:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N addUTR


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
PASA_DIR=pasa_maize

export SINGULARITY_BINDPATH="/path/to/be/mounted"
export SINGULARITY_BIND="/data/folder"

# add UTR using PASA
# http://gaap.hallym.ac.kr/Gaap09

PATH2SQLITE="/scratch/kw68/wt5249/temp/maize/${GENOME}.sqlite"

sed -i -r "s#^(DATABASE=).*#\1$PATH2SQLITE#" /scratch/kw68/wt5249/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/annotCompare.config

# set DATABASE /scratch/kw68/wt5249/sample_data_merged/mysql.confs/alignAssembly.config
# DATABASE=/scratch/kw68/wt5249/temp/mydb_merged.sqlite

cd /scratch/kw68/wt5249/${PASA_DIR}/${GENOME}/sample_data

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/alignAssembly.config -g ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/${GENOME}.genome.fa -P ${SINGULARITY_BIND}/${GENOME}/EVM/EVM.all.gff3

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/annotCompare.config --CPU 2 --TRANSDECODER -A -g ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/${GENOME}.genome.fa -t ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/Trinity.fasta
