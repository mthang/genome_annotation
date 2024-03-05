# Genome Annotation Workflow
1) repeat modeler
2) repeat masking
3) alignment
4) de novo transcriptome assembly
5) genome annotation
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
1)Brachypodium distachyon
2)Oryza sativa subsp. japonica
3)Setaria italica
4)Sorghum bicolor
5)Zea mays


## main steps

