#!/bin/bash

GENOME=Zm-Il14H

CHR_DIR=/maize/${GENOME}/pseudogenes/EVM

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
