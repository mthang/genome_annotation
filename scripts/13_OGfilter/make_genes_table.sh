#!/bin/bash

GROUPS=("core_genes" "shell_genes" "cloud_genes")

for group in core_genes shell_genes cloud_genes;
do
    echo ${group}
    python3 fasta2tab.py ${group}
done
