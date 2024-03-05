#!/bin/bash

#PBS -P your_project_code
#PBS -l ncpus=1
#PBS -l mem=5GB
#PBS -l walltime=1:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N filterPseudo

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

GENOME="genome name/id"
INPUT_DIR=/path/to/pseudogenes/EVM_UTR
OUTPUT_DIR=/path/to/pseudogenes

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr}
    awk '($6 > 0.7 && $11 < 0.0000000001 && $12 > 0.4) {print}' ${INPUT_DIR}/${chr}/output/pgenes/output_pgenes.txt  |  cut -f5,6,11,12,14  |  sort | uniq | grep -v "query" >> ${OUTPUT_DIR}/pseudogene_all_evm_utr_id_frag70_e10_ident40perc.txt
done < ${INPUT_DIR}/chr_scaffold_id.txt


# generate uniq id file
cut -f1 ${OUTPUT_DIR}/pseudogene_all_evm_utr_id_frag70_e10_ident40perc.txt | sort | uniq > ${OUTPUT_DIR}/pseudogene_all_evm_utr_id_frag70_e10_ident40perc_uniq_id.txt

# https://faq.gersteinlab.org/2019/06/03/terms-in-pseudopipe-output-etc/
# frac = fraction of parent gene that matches the pseudogene
# ins = number of insertions in the pseudogene compared to parent sequence
# del = number of deletions in the pseudogene compared to parent sequence
# shift = number of frame shifts in the pseudogene compared to parent sequence
# stop = number of stop codons in the pseudogene compared to parent sequence
# polya = flag indicating the presence or absence of a polyA tail
