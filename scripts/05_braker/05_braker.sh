#!/bin/bash

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=45
#PBS -l mem=80GB
#PBS -l walltime=24:00:0
#PBS -l jobfs=2GB
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N Zm-Il14H_braker

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

### AUGUSTUS BIN DIRECTORY
if [[ -z "${AUGUSTUS_BIN_PATH}" ]];then
       export AUGUSTUS_BIN_PATH=${AUGUSTUS_BIN_PATH}
fi

### AUGUSTUS SCRIPTS DIRECTORY
if [[ -z "${AUGUSTUS_SCRIPTS_PATH}" ]];then
        export AUGUSTUS_SCRIPTS_PATH=${AUGUSTUS_SCRIPTS_PATH}
fi

### BAMTOOLS DIRECTORY
if [[ -z "${BAMTOOLS_PATH}" ]];then
        export BAMTOOLS_PATH=/usr/bin
fi

### GENEMARK DIRECTORY
if [[ -z "${GENEMARK_PATH}" ]];then
        export GENEMARK_PATH=/gmes_linux_64_4
fi

### SAMTOOLS DIRECTORY
if [[ -z "${SAMTOOLS_PATH}" ]];then
        export SAMTOOLS_PATH=/usr/bin
fi

### ALIGNMENT TOOLS DIRECTORY (e.g., GenomeThreader, Spaln, or Exonerate)
if [[ -z "${ALIGNMENT_TOOL_PATH}" ]];then
        export ALIGNMENT_TOOL_PATH=/usr/bin
fi

### PYTHON3 DIRECTORY
if [[ -z "${PYTHON3_PATH}" ]];then
        export PYTHON3_PATH=/usr/bin
fi

### BLAST DIRECTORY
if [[ -z "${BLAST_PATH}" ]];then
        export BLAST_PATH=/usr/bin
fi

export SINGULARITY_BIND="/path/to/be/mounted"

SPECIES=Zm-Il14H
PROTEIN=${SINGULARITY_BINDPATH}/data/master_rename.fasta
ASSEMBLY=${SINGULARITY_BIND}/${SPECIES}/repeatmasker/${SPECIES}.genome.fa.masked
BAM=${SINGULARITY_BIND}/${SPECIES}/hisat


OUTPUT_DIR=/output_dir/${SPECIES}

export PATH=$PATH:/home/564/wt5249/anaconda3/bin

export PROTHINT_PATH=/gmes_linux_64_4/ProtHint/bin

export AUGUSTUS_CONFIG_PATH=/g/data/kw68/analysis/augustus/config
export AUGUSTUS_BIN_PATH=/augustus-3.4.0/bin
export AUGUSTUS_SCRIPTS_PATH=/augustus-3.4.0/scripts

/scratch/kw68/wt5249/braker_new/braker2.sif braker.pl --species=master_Zm-Il14H \
                                           --grass \
                                           --AUGUSTUS_CONFIG_PATH=/g/data/kw68/analysis/augustus/config \
                                           --genome=${ASSEMBLY} \
                                           --prot_seq=${PROTEIN} \
                                           --bam=${BAM}/b73_sorted_all.bam \
                                           --workingdir=${OUTPUT_DIR}/annotation_skipOptimize \
                                           --prg=gth --gth2traingenes --softmasking --gff3 --cores=44 --skipOptimize


