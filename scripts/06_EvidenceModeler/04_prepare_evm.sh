#!/bin/bash

#PBS -P your_project_code
#PBS -l ncpus=2
#PBS -l mem=5GB
#PBS -l walltime=2:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N preEVM

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

export SINGULARITY_BINDPATH="/path/to/be/mounted"
export SINGULARITY_BIND="/data/folder"

OUTPUT_DIR=${SINGULARITY_BINDPATH}/${GENOME}/EVM
mkdir -p ${OUTPUT_DIR}
SING_IMAGE_DIR=/singularity
EVM_SIF=${SING_IMAGE_DIR}/evm-1.1.1.sif

# step 1 - validate GFF3 file format (optional)
#echo "Validating GFF3 file ..."
#${EVM_SIF}  /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl ${PASA_GFF3}
#${EVM_SIF}  /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl ${AUGUSTUS_GFF3}
#${EVM_SIF}  /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl ${GENEMARK_GFF3}
#${EVM_SIF}  /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl ${MERGED_GFF}
#${EVM_SIF}  /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl ${PROTEIN_GFF3}
#echo "End Gff3 File validation."

# ~/BRAKER2/EVM/1.1.1/EvmUtils/partition_EVM_inputs.pl
# step 2 - split contigs

cd ${OUTPUT_DIR}

singularity exec ${EVM_SIF} /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
        --genome ${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked  \
        --gene_predictions ${SINGULARITY_BINDPATH}/${GENOME}/gff/merged/merged.gff \
        --transcript_alignments ${SINGULARITY_BIND}/pasa_maize/${GENOME}/sample_data/${GENOME}.sqlite.pasa_assemblies.gff3 \
        --segmentSize 1000000 --overlapSize 200000 --partition_listing ${SINGULARITY_BINDPATH}/${GENOME}/EVM/partitions_list.out

# step 3 - create weight file
singularity exec ${EVM_SIF} /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/create_weights_file.pl \
        -A ${SINGULARITY_BINDPATH}/${GENOME}/gff/merged/merged.gff \
        -T ${SINGULARITY_BIND}/pasa_maize/${GENOME}/sample_data/${GENOME}.sqlite.pasa_assemblies.gff3 > `pwd`/weights.txt

# step 4 prepare command file
singularity exec ${EVM_SIF}  /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/write_EVM_commands.pl \
        --genome ${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked --weights `pwd`/weights.txt \
        --gene_predictions ${SINGULARITY_BINDPATH}/${GENOME}/gff/merged/merged.gff \
        --transcript_alignments ${SINGULARITY_BIND}/pasa_maize/${GENOME}/sample_data/${GENOME}.sqlite.pasa_assemblies.gff3 \
        --output_file_name evm.out --partitions ${SINGULARITY_BINDPATH}/${GENOME}/EVM/partitions_list.out > ${SINGULARITY_BINDPATH}/${GENOME}/EVM/commands.list

sed -i 's,^,'/g/data/kw68/singularity/evm-1.1.1.sif' ,' commands.list
