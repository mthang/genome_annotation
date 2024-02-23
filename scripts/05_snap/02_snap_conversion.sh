#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=20GB
#PBS -l walltime=5:00:0
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
export SINGULARITY_BIND="/scripts/folder"

SPECIES=Zm-Il14H

ASSEMBLY=${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa
OUTPUT_DIR=${SINGULARITY_BINDPATH}/${SPECIES}/snap

EVM_SIF=/singularity/evm-1.1.1.sif

cd ${SINGULARITY_BINDPATH}/${SPECIES}/snap/

# caution : the snap output might be produced in a pbs log file
#singularity exec /g/data/kw68/singularity/snap-20131129.sif snap ${SINGULARITY_BINDPATH}/${SPECIES}/snap/${SPECIES}.hmm ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa

perl ${SINGULARITY_BIND}/singularity/zff2gff3.pl ${SINGULARITY_BINDPATH}/${SPECIES}/snap/snap.zff > ${SINGULARITY_BINDPATH}/${SPECIES}/snap/snap_zff.gff

singularity exec ${EVM_SIF} /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl ${SINGULARITY_BINDPATH}/${SPECIES}/snap/snap_zff.gff > ${SINGULARITY_BINDPATH}/${SPECIES}/snap/snap.gff
