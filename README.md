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
Kallisto   | [download](https://pachterlab.github.io/kallisto/download)
stringtie  | [download](https://ccb.jhu.edu/software/stringtie/)
bedops     | [download](https://bedops.readthedocs.io/en/latest/index.html)
faSomeRecords | [download](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords)
seqkit | [download](https://bioinf.shenwei.me/seqkit/download/)

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

## Data
### Raw data and Reference genome
The raw data and reference genome download scripts are available [here](https://github.com/mthang/genome_annotation/tree/main/scripts/misc/download)

## Genome Annotation Pipeline
This Genome Annotation pipeline is designed to annotate plant genomes and improve the newly annotated genomes with transcriptome dataset. Genome Annotation pipeline can be summarized in few steps such as repeat detection, gene model prediction, obtain consensus gene model, add utr to the gene model (optional) and retain gene model with gene expression (if transcriptome data is available). Eggnog-mapper is used to perform the functional annotation of the gene model produced in the Genome Annotation Tool, and Orthofinder takes the protein sequences to perform the pangenome analysis.

### Repeat Modeler
#### Input data and Resource
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
#### Input data and Resource
- Reference genome in FASTA format
- Repeats from Repeat modelers in FASTA format
- The Repeat Masking singularity container is used (see the link above)
- PBS script [02_repeat_masker.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/02_repeat_masker/02_repeat_masker.sh) is located in the scripts folder
```
singularity exec ${SINGULARITY_BINDPATH}/tetools_repeat.sif RepeatMasker  -xsmall -pa 24 -gff -rmblast_dir /software/rmblast-2.11.0/bin/ -lib ${SINGULARITY_BINDPATH}/${GENOME}/consensi.fa.classified -dir ${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker ${SINGULARITY_BINDPATH}/${GENOME}/${GENOME}.genome.fa
```

### Aligment
#### Input data and Resource
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
#### Input data and Resource
- Alignment file in BAM format
- The existing genome anontation file (GFF) of the genome of interest
- The installation of stringtie on your computer is required
- PBS script [04_stringtie.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/04_stringtie/04_stringtie.sh) is located in the scripts folder
```
${STRINGTIE} ${HISAT_BAM} -l merged -p 5 -G ${STRINGTIE_GFF_DIR}/${GENOME}/${GENOME}.gff3 -o ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/merged_stringtie.gtf

${STRINGTIE} --merge -p 5 -G ${STRINGTIE_GFF_DIR}/${GENOME}/${GENOME}.gff3 -o ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/final_merged_stringtie.gtf ${STRINGTIE_GFF_DIR}/${GENOME}/stringtie/merged_stringtie.gtf
```

### De novo transcriptome assembly
#### Input data and Resource
- RNAseq data in FASTQ.gz format
- The trinity singularity container is used (see link above)
- PBS script [04_trinity.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/04_trinity/04_trinity.sh) is located in the scripts folder
```
singularity exec ${TRINITY_CONTAINER} Trinity --seqType fq --CPU 50 --max_memory 100G --min_glue 2 --min_kmer_cov 2 --path_reinforcement_distance 75 --group_pairs_distance 250 --min_contig_length 200 --full_cleanup --left ${INPUT_DIR}/${sample}_1.fq.gz --right ${INPUT_DIR}/${sample}_2.fq.gz --output ${FASTQ_DIR}/trinity
```

### Genome annotation 
### BRAKER
#### Input data and Resource
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
#### Input data and Resource
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
#### Input data and Resource
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
#### Input data and Resource
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
#### Input data and Resource
- A list of genome annotation files (gff3) from Braker (i.e augustus and genmark), Fgenesh, SNAP and PASA
- The evidence modeler singularity container is used (see link above)
- PBS scripts are located in [06_EvidenceModeler](https://github.com/mthang/genome_annotation/tree/main/scripts/06_EvidenceModeler) folder.
  - scripts to format input files before running Evidence Modeler
  - [01_convert_gtf.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/01_convert_gtf.sh) is to convert genmark gtf to gff3 format and convert augustus gff3 file to Evidence Modeler gff3 format
  - [02_sort_gff.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/02_sort_gff.sh) is to sort the entry in the gff3 files.
  - [03_merged.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/03_merged.sh) is to merge gff3 from genmark, augustus, Fgenesh and SNAP into a single master gff3 file.
  - [04_prepare_evm.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/04_prepare_evm.sh) is to prepare the executable command to run Evidence Modeler in parallel.
  - [05_run_evm.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/05_run_evm.sh) is used to run evidence modeler program to produce consensus gene models in gff3 format.
  - [06_combine_evm.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/06_EvidenceModeler/06_combine_evm.sh) combines each individual gff3 file into a single gff3 file.

### Gene Structure Annotation and Analysis Using PASA: post Evidence Modeler (add UTR)
#### Input data and Resource
- a final gff3 file from Evidence Modeler
- de novo assembled transcripts from Trinity in FASTA format
- the sqlite file produced in the first PASA run
- a config file (annotCompare.config) from PASA program which can be found in the PASA singularity container
- The PASA singularity container (see link above)
- PBS script [02_run_PASA_utr.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/07_PASA/02_run_PASA_utr.sh) is located in [07_PASA](https://github.com/mthang/genome_annotation/tree/main/scripts/07_PASA) folder.
```
PATH2SQLITE="/scratch/kw68/wt5249/temp/maize/${GENOME}.sqlite"

sed -i -r "s#^(DATABASE=).*#\1$PATH2SQLITE#" /scratch/kw68/wt5249/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/annotCompare.config

cd /scratch/kw68/wt5249/${PASA_DIR}/${GENOME}/sample_data

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/alignAssembly.config -g ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/${GENOME}.genome.fa -P ${SINGULARITY_BIND}/${GENOME}/EVM/EVM.all.gff3

singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/mysql.confs/annotCompare.config --CPU 2 --TRANSDECODER -A -g ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/${GENOME}.genome.fa -t ${SINGULARITY_BINDPATH}/${PASA_DIR}/${GENOME}/sample_data/Trinity.fasta
```

### Pseudogenes detection
### Input data and Resource
- masked reference genome in FASTA format
- final gff3 from Evidence Modeler
- The PseudoPipe tool is used (see link above)
- scripts can be found in [08_pseudogenes](https://github.com/mthang/genome_annotation/tree/main/scripts/08_pseudogenes) folder

### Pre-pseudogenes detection
- Step 1 [01_EVM_GFF2SEQ.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/pre_pseudogenes/01_EVM_GFF2SEQ.sh)
```
GFF3=${SINGULARITY_BINDPATH}/${GENOME}/EVM/EVM.all.gff3

FASTA=${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked

singularity exec ${SING_IMAGE_DIR}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/gff3_file_to_proteins.pl  ${GFF3} ${FASTA} > ${SINGULARITY_BINDPATH}/${GENOME}/EVM/evm_prot.fa
```
- Step 2 [02_make_bed.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/pre_pseudogenes/02_make_bed.sh)
```
DATA_DIR=/genome_maize/${GENOME}/EVM/

OUTPUT_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/exonic_regions

mkdir -p ${OUTPUT_DIR}

gff2bed < ${DATA_DIR}/EVM.all.gff3 > ${OUTPUT_DIR}/evm_exon.bed

awk '{print $1"\t"$4"\t"$2"\t"$3}' ${OUTPUT_DIR}/evm_exon.bed | uniq > ${OUTPUT_DIR}/evm_exon_formatted.bed
```
- Step 3 [03_splitBed.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/pre_pseudogenes/03_splitBed.sh)
```
DATA_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/exonic_regions
cd ${DATA_DIR}

input=${DATA_DIR}/evm_exon_formatted.bed

for chr in `cut -f 1 $input | sort | uniq`;
do
        echo $chr
        grep -w $chr $input > $chr.bed
        ln -s $chr.bed ${chr}_exLocs
done
```
- Step 4 [04_split_chr.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/pre_pseudogenes/04_split_chr.sh)
```
DATA_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/exonic_regions
# GENOME
SPECIES=/genome_maize/${GENOME}/${GENOME}.genome.fa
# Genome from repeat masking
SPECIES_RM=/genome_maize/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked
OUTPUT_DIR=/genome_maize/${GENOME}/pseudogenes/EVM/chr

mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

cut -f1 ${DATA_DIR}/evm_exon.bed | sort | uniq > ${OUTPUT_DIR}/chr_id.txt

FASOMERECORDS=/software/faSomeRecords

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr} > tmp.txt
    ${FASOMERECORDS} ${SPECIES} tmp.txt ${chr}.fa
    ${FASOMERECORDS} ${SPECIES_RM} tmp.txt ${chr}_rm.fa
    rm tmp.txt
done < chr_id.txt
```
- Step 5 [05_create_folder.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/pre_pseudogenes/05_create_folder.sh)
```
DATA_DIR=/genome_maize/${GENOME}/pseudogenes/EVM

CHR_DIR=${DATA_DIR}/chr
EXON_DIR=${DATA_DIR}/exonic_regions

PROTEIN=/genome_maize/${GENOME}/EVM/evm_prot.fa

cd ${DATA_DIR}

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr}
    ls ${CHR_DIR}/${chr}.fa
    ls ${CHR_DIR}/${chr}_rm.fa
    ls ${EXON_DIR}/${chr}.bed
    ls ${PROTEIN}

    mkdir -p ${chr}/input/pep
    mkdir -p ${chr}/input/dna
    mkdir -p ${chr}/input/mysql

        cp ${CHR_DIR}/${chr}.fa ${chr}/input/dna/
    cp ${CHR_DIR}/${chr}_rm.fa ${chr}/input/dna/
    cp ${PROTEIN} ${chr}/input/pep/pep.fa
    cp ${EXON_DIR}/${chr}.bed ${chr}/input/mysql/${chr}_exLocs
done < ${DATA_DIR}/chr/chr_id.txt
```
### Run Pseudogenes detection
### Input and Resource
- Input data will be available after running the scripts in pre-presudogenes detection step
- require python2
- The pseudopipe scripts are not PBS compatible
- Step 1 [01_run_ppipe.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/run_pseudogenes/01_run_ppipe.sh)
```
GENOME="use genome of interest"

## Copy all folders prepared in pre_pseudogenes step to "CHR_DIR" folder below on the compute node
CHR_DIR=/data/maize/${GENOME}/pseudogenes/EVM

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr}
    /software/pgenes/pseudopipe/bin/pseudopipe.sh ${CHR_DIR}/${chr}/output \
    ${CHR_DIR}/${chr}/input/dna/${chr}.fa \
    ${CHR_DIR}/${chr}/input/dna/%s.fa \
    ${CHR_DIR}/${chr}/input/pep/pep.fa \
    ${CHR_DIR}/${chr}/input/mysql/%s_exLocs 0 &
done < ${CHR_DIR}/chr/chr_id.txt
```
- Step 2 [02_part2.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/run_pseudogenes/02_part2.sh)
```
CHR_DIR=/data/maize/${GENOME}/pseudogenes/EVM

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr}
    ls ${CHR_DIR}/${chr}/output/pgenes/minus
    cd ${CHR_DIR}/${chr}/output/pgenes/minus
    source setenvPipelineVars
    python2 /data/software/pgenes/pseudopipe/core/runScripts.py ${chr}
    ls ${CHR_DIR}/${chr}/output/pgenes/plus
    cd ${CHR_DIR}/${chr}/output/pgenes/plus
    source setenvPipelineVars
    python2 /data/software/pgenes/pseudopipe/core/runScripts.py ${chr}

    /data/software/pgenes/pseudopipe/ext/genPgeneResult.sh ${CHR_DIR}/${chr}/output ${CHR_DIR}/${chr}/output/pgenes/output_pgenes.txt
    /data/software/pgenes/pseudopipe/ext/genFullAln.sh ${CHR_DIR}/${chr}/output ${CHR_DIR}/${chr}/output/pgenes/output_pgenes.align.gz
done < ${CHR_DIR}/chr/chr_id.txt
```
### Post-pseudogenes detection
### Input data and Resource
- output_pgenes.txt file from pseudogenes detection
- Step 1 [01_filter_pseudogenes](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/post_pseudogenes/01_filter_pseudogenes.sh) is to filter out the psuedogenes
```
GENOME="genome name or id"
INPUT_DIR=/path/to/pseudogenes/EVM_UTR
OUTPUT_DIR=/path/to/pseudogenes

while IFS="" read -r chr || [ -n "$chr" ]
do
    echo ${chr}
    awk '($6 > 0.7 && $11 < 0.0000000001 && $12 > 0.4) {print}' ${INPUT_DIR}/${chr}/output/pgenes/output_pgenes.txt  |  cut -f5,6,11,12,14  |  sort | uniq | grep -v "query" >> ${OUTPUT_DIR}/pseudogene_all_evm_utr_id_frag70_e10_ident40perc.txt
done < ${INPUT_DIR}/chr_scaffold_id.txt


# generate uniq id file
cut -f1 ${OUTPUT_DIR}/pseudogene_all_evm_utr_id_frag70_e10_ident40perc.txt | sort | uniq > ${OUTPUT_DIR}/pseudogene_all_evm_utr_id_frag70_e10_ident40perc_uniq_id.txt
```
- the output file from Step 1
- the GFF3 with UTR added from PASA
- Step 2 [02_filter_GFF3.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/08_pseudogenes/post_pseudogenes/02_filter_GFF3.sh)
```
# The input of this script is a tab delimiter file (i.e pasa.txt) contains two columns 1) species and 2)PASA num
UPDATE=update

PASA_GFF3_DIR=/genome_maize/pasa

GENOME_DIR=/genome_maize

while IFS="" read -r genome || [ -n "$genome" ]
do
   GENOME=`echo $genome | cut -d" " -f1`
   PASA_NUM=`echo $genome | cut -d" " -f2`
   echo $GENOME
   echo $PASA_NUM

   GFF3=${GENOME}.sqlite.gene_structures_post_PASA_updates.${PASA_NUM}.gff3

   #filter out the DUP genes and keep the PSSD (pseudogenes) only
   grep -v "DUP" ${GENOME_DIR}/${GENOME}/pseudogenes/pseudogene_all_evm_utr_id_frag70_e10_ident40perc.txt > ${GENOME_DIR}/${GENOME}/pseudogenes/pseudogene_all_evm_utr_id_frag70_e10_ident40perc_uniq_id_${UPDATE}.txt

   ID_FILE=${GENOME_DIR}/${GENOME}/pseudogenes/pseudogene_all_evm_utr_id_frag70_e10_ident40perc_uniq_id_${UPDATE}.txt

   mkdir -p ${GENOME_DIR}/${GENOME}/gff/final

   awk '$3=="mRNA" {print}' ${PASA_GFF3_DIR}/${GFF3} | cut -f 9 | sed 's/ID=//g' | sed 's/Parent=//g' | sed 's/;/\t/g' | cut -f 1,2 > ${GENOME_DIR}/${GENOME}/gff/final/mrnaGeneAll.list

   awk 'NR==FNR{a[$0]=$0;next} !($1 in a) {print a[(FNR)],$0}' ${ID_FILE} ${GENOME_DIR}/${GENOME}/gff/final/mrnaGeneAll.list > ${GENOME_DIR}/${GENOME}/gff/final/keep.list
done < pasa.txt
```
### Align RNAseq to CDS and filter out unexpressed CDS
#### Input data and Resource
- Reference genome in FASTA format
- RNASeq data in FASTQ.gz format
- Use PASA singularity container to extract CDS
- Kallisto is used to perform the alignment
- PBS scripts can be found in [09_kallisto](https://github.com/mthang/genome_annotation/tree/main/scripts/09_kallisto) folder.
- Step 1 [01_get_CDS.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/09_kallisto/01_get_CDS.sh)
```
GFF=pasa_pseudogenes_filtered_prekallisto.gff3

mkdir -p ${SINGULARITY_BIND}/${GENOME}/kallisto

# PASA (UTR) after pseudogenes removed
singularity exec ${SINGULARITY_BINDPATH}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/gff3_file_to_proteins.pl ${SINGULARITY_BIND}/${GENOME}/gff/final/${GFF} ${SINGULARITY_BIND}/${GENOME}/${GENOME}.genome.fa CDS > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_wo_pseudogenes.fa
```
- Step 2 [02_run_kallisto.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/09_kallisto/02_run_kallisto.sh)
```
GENOME="genome name or id"

KALLISTO=/software/kallisto/kallisto

INDEX=pasa_wo_pseudogenes.idx

${KALLISTO} index -i ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_wo_pseudogenes.idx ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_wo_pseudogenes.fa

${KALLISTO} quant --threads=20 -i ${SINGULARITY_BIND}/${GENOME}/kallisto/${INDEX} -o ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant ${SINGULARITY_BINDPATH}/b73_merged_1.fq.gz ${SINGULARITY_BINDPATH}/b73_merged_2.fq.gz
```
- Step 3 [03_filter.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/09_kallisto/03_filter.sh)
```
GENOME="genome name or id"

GFF3_DIR=/genome_maize/

ID_FILE=kallisto/pasa_quant/id_no_expression.txt

awk '$5==0 {print}' ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/abundance.tsv > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/cds_no_expression.tsv

awk '$5!=0 {print}' ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/abundance.tsv > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/cds_w_expression.tsv

cut -f1 ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/cds_no_expression.tsv | sort | uniq | grep -v "target_id" > ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/id_no_expression.txt
```
### Statistic of GFF3
#### Input data and Resource
- genome annotation file in GFF3 format
- mikado singularity container is used (see link above)
- PBS scripts [10_mikado](https://github.com/mthang/genome_annotation/tree/main/scripts/10_mikado] is located in the scripts folder
- before filtering out the unexpressed CDS
  - pre-filtering gff3 file run this script [01_mikado_pre_kallisto.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/10_mikado/01_mikado_pre_kallisto.sh)
```
GENOME="genome name or id"
GENOME_DIR=/genome_maize
PASA_NUM=967772

GFF3=${GENOME}.sqlite.gene_structures_post_PASA_updates.${PASA_NUM}.gff3

singularity exec /singularity/mikado-2.3.3.sif mikado util grep ${SINGULARITY_BIND}/${GENOME}/gff/final/keep.list ${SINGULARITY_BINDPATH}/pasa/${GFF3} ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.gff3

singularity exec /singularity/mikado-2.3.3.sif mikado util stats ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.gff3 ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.stats
```
- after filtering out the unexpressed CDS
  - post-filtering gff3 file run this script [02_mikado_post_kallisto.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/10_mikado/02_mikado_post_kallisto.sh)
```
GENOME="genome name or id"

singularity exec /g/data/kw68/singularity/mikado-2.3.3.sif mikado util grep ${SINGULARITY_BIND}/${GENOME}/kallisto/pasa_quant/keep_expression.list ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_prekallisto.gff3 ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_postkallisto.gff3

singularity exec /g/data/kw68/singularity/mikado-2.3.3.sif mikado util stats ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_postkallisto.gff3 ${SINGULARITY_BIND}/${GENOME}/gff/final/pasa_pseudogenes_filtered_postkallisto.stats
```
### Functional Annotation
#### Input data and Resource
- Reference genome in FASTA format
- gff3 file obtained in the post kallisto process
- set up the eggNog database prior to run the functional annotation
- PBS scripts can be found in [11_eggNog](https://github.com/mthang/genome_annotation/tree/main/scripts/11_eggNog]) folder.
- Step 1 [01_get_protein_sequence](https://github.com/mthang/genome_annotation/blob/main/scripts/11_eggNog/01_get_protein_sequence.sh) to obtain protein sequences
```
GENOME="genome name or id"

export SINGULARITY_BINDPATH="/path/to/be/mounted"

SING_IMAGE_DIR=/singularity

GFF3=${SINGULARITY_BINDPATH}/${GENOME}/gff/final/pasa_pseudogenes_filtered_postkallisto.gff3

FASTA=${SINGULARITY_BINDPATH}/${GENOME}/repeatmasker/${GENOME}.genome.fa.masked

singularity exec ${SING_IMAGE_DIR}/pasa_2.5.2.sif /usr/local/src/PASApipeline/misc_utilities/gff3_file_to_proteins.pl  ${GFF3} ${FASTA} > ${SINGULARITY_BINDPATH}/${GENOME}/gff/final/final_prot.fa
```
- Step 2 [02_run_eggnog.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/11_eggNog/02_run_eggnog.sh) annotate protein sequences
```
GENOME="genome name or id"

SCRATCH_DIR=/scratch
TEMP_DIR=/tmp

singularity exec ${SINGULARITY_BIND}/singularity/eggnog-mapper_2.1.9.sif emapper.py --data_dir /g/data/kw68/data/eggnogDB -m 'diamond' -i ${SINGULARITY_BIND}/genome_maize/${GENOME}/gff/final/final_prot.fa --itype 'proteins' --matrix 'BLOSUM62' --gapopen 11 --gapextend 1 --sensmode sensitive --dmnd_iterate no --score 0.001  --seed_ortholog_evalue 0.001 --target_orthologs=all --go_evidence=non-electronic --no_file_comments --report_orthologs  --output=${SINGULARITY_BIND}/genome_maize/${GENOME}/gff/final/eggNog --cpu 12 --scratch_dir ${SCRATCH_DIR} --temp_dir ${TEMP_DIR}
```
### Pangenome Analysis
#### Input data and Resource
- eggNog Functional annotation output file
- protein sequences fasta file used in eggNog step
- Orthofinder singularity container is used (see link above)
- scripts can be found in [12_orthofinder]https://github.com/mthang/genome_annotation/tree/main/scripts/12_orthofinder) folder.
- Step 1 [01_reformat_fasta_header.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/12_orthofinder/01_reformat_fasta_header.sh)
```
GENOME_DIR=/genome_maize

while IFS="" read -r genome || [ -n "$genome" ]
do
   GENOME=`echo $genome | cut -d" " -f1`
   echo $GENOME
   cd ${GENOME_DIR}/${GENOME}/gff/final

   cut -f1-2,8-9 eggNog.emapper.annotations | grep -v "^#" > header_format.txt
   awk -v var="$GENOME" 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"|"var"|"$1"|"$3"|"$4}' header_format.txt > header_format_tab.txt

   sed '/^>/ s/ .*//' final_prot.fa > final_prot_update.fa
   cut -f1 header_format_tab.txt > fasta_id.txt
   faSomeRecords final_prot_update.fa fasta_id.txt ${GENOME}_prot.fa

   # sed -i '/^>/ s/ .*//' test.fa
   #/g/data/kw68/software/seqkit replace -p "(.*)" --replacement "{kv}" --kv-file tmp_header.txt test.fa

   /g/data/kw68/software/seqkit replace -p "(.+)" -r '{kv}' -k header_format_tab.txt ${GENOME}_prot.fa > ${GENOME}_filtered.fa

   sed -i 's/\ /\_/g' ${GENOME}_filtered.fa

done < pasa.txt
```
- a list of protein sequences fasta file 
- Step 2 [02_run_orthofinder_all.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/12_orthofinder/02_run_orthofinder_all.sh)
```
# make a symbolic link to the *_filtered.fa (protein sequences)
DATA_DIR=${SINGULARITY_BIND}/data/orthofinder_all

#memory exceeded issue when the number species is growing
singularity exec ${SINGULARITY_BIND}/singularity/OrthoFinder-2.5.4.sif orthofinder -t 30 -f ${DATA_DIR}
```

### Orthologues Genes
#### Input data and Resource
- Orthogroups.GeneCount.tsv file from Orthofinder
- Orthogroup Sequences in Orthogroups_Sequences folder
- Custom scripts in [13_GOfilter](https://github.com/mthang/genome_annotation/tree/main/scripts/13_OGfilter) folder is used
- make changes in [runOGfilter.sh](https://github.com/mthang/genome_annotation/blob/main/scripts/13_OGfilter/runOGfilter.sh) to subset core, shell and cloud genes.
  - Core genes - a unique gene found in all species (N);
     - set MIN_SPECIES = 0.92 in runOGfilter.sh
     - use [OGFilter_max.py](https://github.com/mthang/genome_annotation/blob/main/scripts/13_OGfilter/OGFilter_max.py)
  - Shell genes - a unique gene found in all species (N - 1);
     - set MIN_SPECIES = 0.15 in runOGfilter.sh
     - use [OGFilter_min_max.py](https://github.com/mthang/genome_annotation/blob/main/scripts/13_OGfilter/OGFilter_min_max.py)  
  - Cloud genes - a unique gene found in a single species;
     - set MIN_SPECIES = 0.15 in runOGfilter.sh
     - use [OGFilter_min.py](https://github.com/mthang/genome_annotation/blob/main/scripts/13_OGfilter/OGFilter_min.py)
- this [script](https://github.com/mthang/genome_annotation/blob/main/scripts/13_OGfilter/make_genes_table.sh) is to convert the output file produced by runOGfilter.sh to a tabular file with the following columns:
   - Ortholog_Filename
   - Taxonomy
   - Gene
   - Sequence_Name
   - Description
   - Preferred_Name
   - Protein_Sequence
## Reference
