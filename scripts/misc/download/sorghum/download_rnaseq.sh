#!/bin/bash


#URL to the ftp server
URL=https://ftp.cngb.org/pub/CNSA/data3/CNP0001440

#RNAseq data
FILES=(CNS0314252/CNX0351653/CNR0431540/353_1.fq.gz 
CNS0314254/CNX0351654/CNR0431541/IS3614-3_1.fq.gz 
CNS0314255/CNX0351656/CNR0431543/IS929_1.fq.gz 
CNS0314258/CNX0351655/CNR0431542/IS8525_1.fq.gz 
CNS0314259/CNX0351657/CNR0431544/Ji2731_1.fq.gz 
CNS0314260/CNX0351658/CNR0431545/PI525695_1.fq.gz)

# reverse_end=`echo ${FASTQ_DIR}/${p} | sed 's/_1/_2/g'`

for f in ${FILES[@]}
do
   if [[ "${f}" != *.fa.gz && "${f}" != *.gfa && "${f}" != *.txt ]];then
       echo ${URL}/${f}
       wget -nH --cut-dirs=3 -x ${URL}/${f}
       if [[ "${f}" == *_1.fq.gz ]]; then
          reverse_end=`echo ${f} | sed 's/_1/_2/g'`
          echo ${URL}/${reverse_end}
          wget -nH --cut-dirs=3 -x ${URL}/${reverse_end}
       fi
   else
       echo ${URL}/${f}
       echo "downloading misc"
       wget -nH --cut-dirs=3 -x ${URL}/${f}
   fi
done
