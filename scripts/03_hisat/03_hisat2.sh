#!/bin/bash


##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

#PBS -P your_project_code
#PBS -l ncpus=20
#PBS -l mem=100GB
#PBS -l walltime=20:00:0
#PBS -l storage=scratch/kw68+gdata/kw68
#PBS -N hi_Zm-Il14H

module load singularity

module load samtools/1.12

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

export SINGULARITY_BINDPATH="/path/to/be/mounted/"
export SINGULARITY_BIND="/data/path/to/be/mounted"

cpu=20
GENOME=Zm-Il14H
#BAM file
SAMPLE=(b73_merged)

for sample in ${SAMPLE[@]}
do
sampleName=`echo $sample | sed 's/\_merged*.//g'`

singularity exec ${SINGULARITY_BINDPATH}/hisat2.sif hisat2 -p ${cpu} -x ${SINGULARITY_BIND}/${GENOME}/${GENOME}.genome.fa \
       -1 ${SINGULARITY_BINDPATH}/${sample}_1.fq.gz \
       -2 ${SINGULARITY_BINDPATH}/${sample}_2.fq.gz \
       -S ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.sam \
       --novel-splicesite-outfile ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.junctions \
       --rna-strandness RF \
       --dta \
       -t &> ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.log

samtools view \
       --threads ${cpu} \
       -b \
       -o ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.bam \
       ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.sam

# multi thread command when running on a single machine 
# samtools view -b -F 4 -@ 44 ${SINGULARITY_BIND}/analysis_merged/hisat/${sample}.bam > ${SINGULARITY_BIND}/analysis_merged/hisat/${sample}_mapped.bam


# Resource Allocation = 4G (memory) x 20 (CPU) = 80G (memory) <- do not exceed 100G (memory), so the maximum number of cpu is 25 (4x25 =100G)

samtools sort \
      -m 4G \
      -o ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}_sorted_all.bam \
      -T ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}_temp \
      --threads ${cpu} \
      ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.bam

      rm ${SINGULARITY_BIND}/${GENOME}/hisat/${sampleName}.sam
done
