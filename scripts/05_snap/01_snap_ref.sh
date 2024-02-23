#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=30GB
#PBS -l walltime=30:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Zm-Il14H_snap
#PBS -e error.txt
#PBS -o output.txt

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
export SINGULARITY_BIND="/tool/folder"

SPECIES=Zm-Il14H
ASSEMBLY=${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa

OUTPUT_DIR=${SINGULARITY_BINDPATH}/${SPECIES}/snap

# step 1
perl ${SINGULARITY_BIND}/singularity/gff3_to_zff.pl ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.gff3 > ${OUTPUT_DIR}/${SPECIES}.zff

cd ${OUTPUT_DIR}

# step 2
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom -validate ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  > ${OUTPUT_DIR}/${SPECIES}.validate

singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa -gene-stats > ${OUTPUT_DIR}/gene-stats.log 2>&1

singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  -validate > ${OUTPUT_DIR}/validate.log 2>&1
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  -categorize 1000 > ${OUTPUT_DIR}/categorize.log 2>&1
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/uni.ann ${OUTPUT_DIR}/uni.dna -export 1000 -plus > ${OUTPUT_DIR}/uni-plus.log 2>&1


mkdir params
cd params

singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif forge ${OUTPUT_DIR}/export.ann ${OUTPUT_DIR}/export.dna > ${OUTPUT_DIR}/forge.log 2>&1

perl ${SINGULARITY_BIND}/singularity/hmm-assembler.pl ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  ${OUTPUT_DIR}/params/ > ${OUTPUT_DIR}/${SPECIES}.hmm

# caution : the snap output might be produced in a pbs log file
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif snap -gff ${OUTPUT_DIR}/${SPECIES}.hmm ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa > ${OUTPUT_DIR}/snap.zff
