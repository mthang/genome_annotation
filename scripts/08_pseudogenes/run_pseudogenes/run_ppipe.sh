
#!/bin/bash

GENOME=Zm-Il14H

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
