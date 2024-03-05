#!/bin/bash

## uncomment GENOME and PASA_NUM for individual species
#GENOME=9311
#PASA_NUM=3491083

# The input of this script is a tab delimiter file (i.e pasa.txt) contains two columns 1) species and 2)PASA num
UPDATE=update

#PASA_GFF3_DIR=/genome_maize/pasa/${GENOME}
PASA_GFF3_DIR=/genome_maize/pasa
#cd ${PASA_GFF3_DIR}

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
