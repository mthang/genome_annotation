#!/bin/bash

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
   /home/564/wt5249/software/faSomeRecords final_prot_update.fa fasta_id.txt ${GENOME}_prot.fa

   # sed -i '/^>/ s/ .*//' test.fa
   #/g/data/kw68/software/seqkit replace -p "(.*)" --replacement "{kv}" --kv-file tmp_header.txt test.fa

   /g/data/kw68/software/seqkit replace -p "(.+)" -r '{kv}' -k header_format_tab.txt ${GENOME}_prot.fa > ${GENOME}_filtered.fa

   sed -i 's/\ /\_/g' ${GENOME}_filtered.fa

done < pasa.txt

# run bash parser.sh > newtest.fa
# check empty annotations grep "> *$" newtest.fa | wc -l
