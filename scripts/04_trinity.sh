#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --job-name=Trinity
#SBATCH --time=96:00:00
#SBATCH --partition=general
#SBATCH --account=your_account
#SBATCH -o slurm_trinity.output
#SBATCH -e slurm_trinity.error

###########
# Due to the limited walltime (48 hours) on NCI HPC, this SLURM script is used to run on BUNYA machine for de novo assembly
##########
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=C

FASTQ_DIR=/merged
INPUT_DIR=/merged

TRINITY_CONTAINER=/software/trinityrnaseq.v2.13.2.simg

sample=b73_merged

export _JAVA_OPTIONS="-Xmx20g"

singularity exec ${TRINITY_CONTAINER} Trinity --seqType fq --CPU 50 --max_memory 100G --min_glue 2 --min_kmer_cov 2 --path_reinforcement_distance 75 --group_pairs_distance 250 --min_contig_length 200 --full_cleanup --left ${INPUT_DIR}/${sample}_1.fq.gz --right ${INPUT_DIR}/${sample}_2.fq.gz --output ${FASTQ_DIR}/trinity
