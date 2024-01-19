#Genome Annotation Tools

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
14) https://github.com/PASApipeline/PASApipeline/blob/master/docs/index.asciidoc

## Singularity Image
Tool  | Singularity Container | URL
-----------------|--------------------|------------
Repeat Modeler   | tetools_repeat.sif | [Github](https://github.com/Dfam-consortium/TETools)
Repeat Masker    | tetools_repeat.sif | [Github](https://github.com/Dfam-consortium/TETools)
Mikado           | mikado-2.3.3.sif   | [Link](https://quay.io/repository/biocontainers/mikado?tab=history)
Snap             | snap-20131129.sif  |
braker2          | braker2.sif        | [Link](https://quay.io/repository/biocontainers/braker2?tab=tags)
hisat2           | hisat2.sif         | [Link](https://quay.io/repository/biocontainers/hisat2?tab=tags)
fgenesh          | fgenesh.sif        | private
Evidence Modeler | evm-1.1.1.sif      | [Link](https://quay.io/repository/biocontainers/evidencemodeler?tab=tags)
PASA             | pasa_2.5.2.sif     | [Link](https://quay.io/repository/biocontainers/pasa?tab=tags)
Eggnog-mapper    | eggnog-mapper_2.1.9.sif | [Link](https://quay.io/repository/biocontainers/eggnog-mapper?tab=tags)
OrthoFinder      | OrthoFinder-2.5.4.sif   | [Link](https://quay.io/repository/biocontainers/orthofinder?tab=tags)

