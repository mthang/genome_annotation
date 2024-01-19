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
---------------| -------------
Repeat Modeler | tetools_repeat.sif
Repeat Masker  | tetools_repeat.sif
Mikado         | mikado-2.3.3.sif

* snap-20131129.sif <br/>
* evm-1.1.1.sif 
* hisat2.sif
* pasa_2.5.2.sif
* fgenesh.sif
* braker2.sif
* eggnog-mapper_2.1.9.sif
* OrthoFinder-2.5.4.sif
