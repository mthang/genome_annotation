#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=20
#PBS -l mem=20GB
#PBS -l walltime=30:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N pasa

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

mkdir ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif cp -r /usr/local/src/PASApipeline/sample_data/ ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/

# set DATABASE /scratch/kw68/wt5249/sample_data_merged/mysql.confs/alignAssembly.config
# DATABASE=/scratch/kw68/wt5249/temp/mydb_merged.sqlite

PATH2SQLITE="/scratch/kw68/wt5249/temp/maize/${GENOME}.sqlite"

sed -i -r "s#^(DATABASE=).*#\1$PATH2SQLITE#" ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/mysql.confs/alignAssembly.config

# Merged RNASeq

cp ${SINGULARITY_BIND}/data/trinity/maize/trinity.Trinity.fasta ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/Trinity.fasta
cp ${SINGULARITY_BIND}/genome_maize/${GENOME}/${GENOME}.genome.fa ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/

# Change directory
cd ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/accession_extractor.pl < ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/Trinity.fasta > ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/tdn.accs

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/mysql.confs/alignAssembly.config --trans_gtf ${SINGULARITY_BIND}/genome_maize/${GENOME}/stringtie/merged_stringtie.gtf -C -R --CPU 40 --ALIGNER gmap -g ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/${GENOME}.genome.fa -t ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/Trinity.fasta --TDN ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/tdn.accs
