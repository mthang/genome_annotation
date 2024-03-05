# Genome Annotation Workflow
1) Repeat Modeler
2) Repeat Masking
3) Alignment
4) De novo transcriptome assembly
5) Genome annotation
   1. Braker
   2. Snap
   3. Fgenesh
7) Gene Structure Annotation and Analysis Using PASA: pre Evidence Modeler
8) Evidence Modeler
9) Gene Structure Annotation and Analysis Using PASA: post Evidence Modeler (add UTR)
10) Pseudogenes detection
11) Expression of CDS
12) Filter unexpressed CDS from GFF3
13) GFF3 statistics
14) Functional Annotation
15) Pangenome analysis
16) https://github.com/PASApipeline/PASApipeline/blob/master/docs/index.asciidoc

## Singularity Image
Tool  | Singularity Container | URL
-----------------|--------------------|------------
Repeat Modeler   | tetools_repeat.sif | [Github](https://github.com/Dfam-consortium/TETools)
Repeat Masker    | tetools_repeat.sif | [Github](https://github.com/Dfam-consortium/TETools)
Mikado           | mikado-2.3.3.sif   | [Link](https://quay.io/repository/biocontainers/mikado?tab=history)
Trinity          | trinityrnaseq.v2.13.2.simg | [Link](https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/archived/trinityrnaseq.v2.13.2.simg)
Snap             | snap-20131129.sif  | [Link](https://hub.docker.com/r/biocontainers/snap/tags)
braker2          | braker2.sif        | [Link](https://quay.io/repository/biocontainers/braker2?tab=tags)
hisat2           | hisat2.sif         | [Link](https://quay.io/repository/biocontainers/hisat2?tab=tags)
fgenesh          | fgenesh.sif        | private
Evidence Modeler | evm-1.1.1.sif      | [Link](https://quay.io/repository/biocontainers/evidencemodeler?tab=tags)
PASA             | pasa_2.5.2.sif     | [Link](https://quay.io/repository/biocontainers/pasa?tab=tags)
Eggnog-mapper    | eggnog-mapper_2.1.9.sif | [Link](https://quay.io/repository/biocontainers/eggnog-mapper?tab=tags)
OrthoFinder      | OrthoFinder-2.5.4.sif   | [Link](https://quay.io/repository/biocontainers/orthofinder?tab=tags)

## Binary Tool
Tool | URL
-----------|-----------
Kallisto   | [Link](https://pachterlab.github.io/kallisto/download)
stringtie  | [Link](https://ccb.jhu.edu/software/stringtie/)

## Single Threaded Tool
Tool | URL
-----------|-----------
PseudoPipe | [Link](http://pseudogene.org/pseudopipe/)

## Homologue Database
Five non-redundant databases are downloaded from [Uniprot](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/) 
- Brachypodium distachyon-
- Oryza sativa subsp. japonica
- Setaria italica
- Sorghum bicolor
- Zea mays

## Genome Annotation Pipeline
This Genome Annotation pipeline is designed to annotate plant genomes and improve the newly annotated genomes with transcriptome dataset. Genome Annotation pipeline can be summarized in few steps such as repeat detection, gene model prediction, obtain consensus gene model, add utr to the gene model (optional) and retain gene model with gene expression (if transcriptome data is available). Eggnog-mapper is used to perform the functional annotation of the gene model produced in the Genome Annotation Tool, and Orthofinder takes the protein sequences to perform the pangenome analysis.

### Repeat Modeler
#### Raw data and Resource
- Reference genome / de novo assembled genome in FASTA format
- The Repeat Modeler singularity container is used (see the link above)
- PBS script [01_repeatmodeler.sh](https://github.com/mthang/genome_annotation/tree/main/scripts/01_repeat_modeler) is located in the scripts folder
```
# Step 1 - index reference genome fasta file 
singularity exec ${SINGULARITY_BINDPATH}/tetools_repeat.sif BuildDatabase -name ${SPECIES} -engine ncbi ${SPECIES}.fa

# Step 2 - build the repeat model (putative repeats) of the input reference genome
singularity exec ${SINGULARITY_BINDPATH}/tetools_repeat.sif RepeatModeler -engine ncbi -threads 24 -database ${SPECIES}

Where ${SINGULARITY_BINDPATH} is the variable defined the location of tool folder 
```

### Repeat Masker
#### Raw data and Resource
- Reference genome in FASTA format
- Repeats from Repeat modelers in FASTA format
- The Repeat Masking singularity container is used (see the link above)
- PBS script [02_repeat_masker.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/02_repeat_masker/02_repeat_masker.sh) is located in the scripts folder
```
singularity exec ${SINGULARITY_BINDPATH}/tetools_repeat.sif RepeatMasker  -xsmall -pa 24 -gff -rmblast_dir /software/rmblast-2.11.0/bin/ -lib ${SINGULARITY_BINDPATH}/${GENOME}/consensi.fa.classified -dir ${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker ${SINGULARITY_BINDPATH}/${GENOME}/${GENOME}.genome.fa
```

### Aligment
- Reference genome in FASTA format
- RNAseq data (i.e single- or paired-end) in FASTQ.gz format
- Hisat2 and samtools for alignment and sorting BAM file respectively.
- PBS script [03_hisat2.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/03_hisat/03_hisat2.sh) is located in the scripts folder
```
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
```

## Reference
