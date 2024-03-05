# Genome Annotation Workflow
1) Repeat Modeler
2) Repeat Masking
3) Alignment
4) Merge BAM with GFF file
5) De novo transcriptome assembly
6) Genome annotation
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
#### Raw data and Resource
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
###  Merge BAM with GFF file
#### Raw data and Resource
- Alignment file in BAM format
- The existing genome anontation file (GFF) of the genome of interest
- The installation of stringtie on your computer is required
- PBS script [04_stringtie.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/04_stringtie/04_stringtie.sh) is located in the scripts folder
```
${STRINGTIE} ${HISAT_BAM} -l merged -p 5 -G ${STRINGTIE_GFF_DIR}/${GENOME}/${GENOME}.gff3 -o ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/merged_stringtie.gtf

${STRINGTIE} --merge -p 5 -G ${STRINGTIE_GFF_DIR}/${GENOME}/${GENOME}.gff3 -o ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/final_merged_stringtie.gtf ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/merged_stringtie.gtf
```

### De novo transcriptome assembly
#### Raw data and Resource
- RNAseq data in FASTQ.gz format
- The trinity singularity container is used (see link above)
- PBS script [04_trinity.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/04_trinity/04_trinity.sh) is located in the scripts folder
```
singularity exec ${TRINITY_CONTAINER} Trinity --seqType fq --CPU 50 --max_memory 100G --min_glue 2 --min_kmer_cov 2 --path_reinforcement_distance 75 --group_pairs_distance 250 --min_contig_length 200 --full_cleanup --left ${INPUT_DIR}/${sample}_1.fq.gz --right ${INPUT_DIR}/${sample}_2.fq.gz --output ${FASTQ_DIR}/trinity
```

### Genome annotation 
### BRAKER
#### Raw data and Resource
- Reference genome in FASTA format
- RNAseq data in BAM format
- Homologues Sequences in FASTA format
- The braker singularity container is used (see link above)
- PBS script [05_braker.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/05_braker/05_braker.sh) is located in the scripts folder
```
braker2.sif braker.pl --species=master_Zm-Il14H \
                                           --grass \
                                           --AUGUSTUS_CONFIG_PATH=/g/data/kw68/analysis/augustus/config \
                                           --genome=${ASSEMBLY} \
                                           --prot_seq=${PROTEIN} \
                                           --bam=${BAM}/b73_sorted_all.bam \
                                           --workingdir=${OUTPUT_DIR}/annotation_skipOptimize \
                                           --prg=gth --gth2traingenes --softmasking --gff3 --cores=44 --skipOptimize
```

### Fgenesh
#### Raw data and Resource
- Repeat masked reference genome in FASTA format
- a fgenesh [config file](https://github.com/mthang/genome_annotation/blob/main/scripts/05_fgenesh/plant.cfg) which defines path to the reference files (i.e GENE_PARAM, PIPE_PARAM, PROTEIN_DB and BLASTP) and tools.
   - GENE_PARAM - gene prediction parameters (i.e gene matrix provided by fgenesh)
   - PIPE_PARAM - location of fgenesh parameters files 
   - PROTEIN_DB - protein database
   - BLASTP - path to executable blastp program
- The fgenesh singulartiy image 
- PBS script [05_fgenesh.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/05_fgenesh/05_fgenesh.sh) is located in the scripts folder
```
cpu=20
GENOME=GenomeOfInterest

cd ${SINGULARITY_BIND}/${GENOME}/fgenesh
mkdir gff3
mkdir results
mkdir contigs

#1)  make sequence list
singularity exec ${SINGULARITY_BINDPATH}/singularity/fgenesh.sif split_multi_fasta.pl ${SINGULARITY_BIND}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked \
       -name seq_id \
       -dir ${SINGULARITY_BIND}/${GENOME}/fgenesh/contigs \
       -ext fa \
       -col 60 \
       -mklist ${SINGULARITY_BIND}/${GENOME}/fgenesh/contigs.list

#2) make sequence list with softmasked sequence
singularity exec  ${SINGULARITY_BINDPATH}/singularity/fgenesh.sif split_multi_fasta.pl ${SINGULARITY_BIND}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked \
        -name seq_id \
        -dir ${SINGULARITY_BIND}/${GENOME}/fgenesh/contigs \
        -ext fa.masked \
        -col 60 \
        -mklist ${SINGULARITY_BIND}/${GENOME}/fgenesh/contigsN.list

#3) run fgenesh
singularity exec  ${SINGULARITY_BINDPATH}singularity/fgenesh.sif run_pipe.pl ${SINGULARITY_BIND}/${GENOME}/fgenesh/plant.cfg \
        -l ${SINGULARITY_BIND}/${GENOME}/fgenesh/contigs/contigs.list \
        -m ${SINGULARITY_BIND}/${GENOME}/fgenesh/contigs/contigsN.list \
        -d ${SINGULARITY_BIND}/${GENOME}/fgenesh/results \
        -w ${SINGULARITY_BIND}/${GENOME}/fgenesh/tmp 2> run_plant.log

#4 run convert resn3 to gff3 format
singularity exec  ${SINGULARITY_BINDPATH}/singularity/fgenesh.sif run_fgenesh_2_gff3.pl \
        ${SINGULARITY_BIND}/${GENOME}/fgenesh/results \
        ${SINGULARITY_BIND}/${GENOME}/fgenesh/gff3 -sort -print_exons
```

### SNAP
#### Raw data and Resource
- The existing genome annnotation file (gff3) of the genome of interest downloaded from the public repository (i.e NCBI)
- The snap singulariy container (see link above)
- PBS script is located in the scripts folder
   - convert reference to snap format and run snap [01_snap_ref.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/05_snap/01_snap_ref.sh)
   - convert snap output file (zff) to gff file format [02_snap_conversionh](https://github.com/mthang/genome_annotation/blob/main/scripts/05_snap/02_snap_conversion.sh). This has been implemented in the script 01_snap_ref.sh
```
ASSEMBLY=${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa

OUTPUT_DIR=${SINGULARITY_BINDPATH}/${SPECIES}/snap

# step 1 -transform the existing genome annotation file (gff3) to the snap format (zff)
perl ${SINGULARITY_BIND}/singularity/gff3_to_zff.pl ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.gff3 > ${OUTPUT_DIR}/${SPECIES}.zff

cd ${OUTPUT_DIR}

# step 2 - validate the zff file
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom -validate ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  > ${OUTPUT_DIR}/${SPECIES}.validate

singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa -gene-stats > ${OUTPUT_DIR}/gene-stats.log 2>&1

singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  -validate > ${OUTPUT_DIR}/validate.log 2>&1
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/${SPECIES}.zff ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  -categorize 1000 > ${OUTPUT_DIR}/categorize.log 2>&1
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif fathom ${OUTPUT_DIR}/uni.ann ${OUTPUT_DIR}/uni.dna -export 1000 -plus > ${OUTPUT_DIR}/uni-plus.log 2>&1

# step 3 - create paramater files
mkdir params
cd params

singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif forge ${OUTPUT_DIR}/export.ann ${OUTPUT_DIR}/export.dna > ${OUTPUT_DIR}/forge.log 2>&1

perl ${SINGULARITY_BIND}/singularity/hmm-assembler.pl ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa  ${OUTPUT_DIR}/params/ > ${OUTPUT_DIR}/${SPECIES}.hmm

# step 4 - run snap
# caution : the snap output might be produced in a pbs log file
singularity exec ${SINGULARITY_BIND}/singularity/snap-20131129.sif snap -gff ${OUTPUT_DIR}/${SPECIES}.hmm ${SINGULARITY_BINDPATH}/${SPECIES}/${SPECIES}.genome.fa > ${OUTPUT_DIR}/snap.zff
```

### Gene Structure Annotation and Analysis Using PASA: pre Evidence Modeler
#### Raw data and Resource
- de novo assembled transcripts from trinity in FASTA format
- a config file (alignAssembly.config) from PASA program which can be found in the PASA singularity container
- The PASA singularity container (see link above)
- PBS script [01_run_PASA.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/07_PASA/01_run_PASA.sh) is located in the scripts folder
```
GENOME=Zm-Il14H

mkdir ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif cp -r /usr/local/src/PASApipeline/sample_data/ ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/

PATH2SQLITE="/scratch/kw68/wt5249/temp/maize/${GENOME}.sqlite"

# replace the default path defined in the alignAssembly.config with PATH2SQLITE
sed -i -r "s#^(DATABASE=).*#\1$PATH2SQLITE#" ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/mysql.confs/alignAssembly.config

# Merged RNASeq

cp ${SINGULARITY_BIND}/data/trinity/maize/trinity.Trinity.fasta ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/Trinity.fasta
cp ${SINGULARITY_BIND}/genome_maize/${GENOME}/${GENOME}.genome.fa ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/

# Change directory
cd ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/accession_extractor.pl < ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/Trinity.fasta > ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/tdn.accs

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/mysql.confs/alignAssembly.config --trans_gtf ${SINGULARITY_BIND}/genome_maize/${GENOME}/stringtie/merged_stringtie.gtf -C -R --CPU 40 --ALIGNER gmap -g ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/${GENOME}.genome.fa -t ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/Trinity.fasta --TDN ${SINGULARITY_BINDPATH}/pasa_maize/${GENOME}/sample_data/tdn.accs
```

###  Evidence Modeler
#### Raw data and Resource
- A list of genome annotation files (gff3) from Braker (i.e augustus and genmark), Fgenesh, SNAP and PASA
- The evidence modeler singularity container is used (see link above)
- PBS scripts [06_EvidenceModeler](https://github.com/mthang/genome_annotation/tree/main/scripts/06_EvidenceModeler) is located in the scripts folder
  - scripts to format input files before running Evidence Modeler
  - [01_convert_gtf.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/01_convert_gtf.sh) is to convert genmark gtf to gff3 format and convert augustus gff3 file to Evidence Modeler gff3 format
  - [02_sort_gff.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/02_sort_gff.sh) is to sort the entry in the gff3 files.
  - [03_merged.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/03_merged.sh) is to merge gff3 from genmark, augustus, Fgenesh and SNAP into a single master gff3 file.
  - [04_prepare_evm.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/04_prepare_evm.sh) is to prepare the executable command to run Evidence Modeler in parallel.
  - [05_run_evm.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/05_run_evm.sh) is used to run evidence modeler program to produce consensus gene models in gff3 format.
  - [06_combine_evm.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/06_combine_evm.sh) combines each individual gff3 file into a single gff3 file.

### Gene Structure Annotation and Analysis Using PASA: post Evidence Modeler (add UTR)
#### Raw data and Resource
- a final gff3 file from Evidence Modeler
- de novo assembled transcripts from Trinity in FASTA format
- the sqlite file produced in the first PASA run
- a config file (annotCompare.config) from PASA program which can be found in the PASA singularity container
- The PASA singularity container (see link above)
- PBS script [02_run_PASA_utr.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/07_PASA/02_run_PASA_utr.sh) is located in the scripts folder
```
PATH2SQLITE="/scratch/kw68/wt5249/temp/maize/${GENOME}.sqlite"

sed -i -r "s#^(DATABASE=).*#\1$PATH2SQLITE#" /scratch/kw68/wt5249/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/annotCompare.config

cd /scratch/kw68/wt5249/${PASA_DIR}/${GENOME}/sample_data

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/alignAssembly.config -g ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/${GENOME}.genome.fa -P ${SINGULARITY_BIND}/${GENOME}/EVM/EVM.all.gff3

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/annotCompare.config --CPU 2 --TRANSDECODER -A -g ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/${GENOME}.genome.fa -t ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/Trinity.fasta
```

## Reference
