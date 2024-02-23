#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -l ncpus=12
#PBS -l mem=30GB
#PBS -l walltime=2:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Zm-Il14H
#PBS -e eggNog_error.txt
#PBS -o eggNog_output.txt

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

export SINGULARITY_BIND="/path/to/be/mounted"
export EGGNOG_DATA_DIR=/data/eggnogDB

GENOME=Zm-Il14H

SCRATCH_DIR=/scratch
TEMP_DIR=/tmp

singularity exec ${SINGULARITY_BIND}/singularity/eggnog-mapper_2.1.9.sif emapper.py --data_dir /g/data/kw68/data/eggnogDB -m 'diamond' -i ${SINGULARITY_BIND}/genome_maize/${GENOME}/gff/final/final_prot.fa --itype 'proteins' --matrix 'BLOSUM62' --gapopen 11 --gapextend 1 --sensmode sensitive --dmnd_iterate no --score 0.001  --seed_ortholog_evalue 0.001 --target_orthologs=all --go_evidence=non-electronic --no_file_comments --report_orthologs  --output=${SINGULARITY_BIND}/genome_maize/${GENOME}/gff/final/eggNog --cpu 12 --scratch_dir ${SCRATCH_DIR} --temp_dir ${TEMP_DIR}
