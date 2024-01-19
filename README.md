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
Tool  | Singularity Container 
-----------------|--------------
Repeat Modeler   | tetools_repeat.sif
Repeat Masker    | tetools_repeat.sif
Mikado           | mikado-2.3.3.sif
Snap             | snap-20131129.sif
braker2          | braker2.sif
hisat2           | hisat2.sif
fgenesh          | fgenesh.sif
Evidence Modeler | evm-1.1.1.sif
PASA             | pasa_2.5.2.sif
Eggnog-mapper    | eggnog-mapper_2.1.9.sif
OrthoFinder      | OrthoFinder-2.5.4.sif

