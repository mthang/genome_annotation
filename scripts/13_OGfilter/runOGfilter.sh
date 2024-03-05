#!/bin/bash

GENECOUNT=Results_Jan03/Orthogroups/Orthogroups.GeneCount.tsv
ORTHO_SEQ=Results_Jan03/Orthogroup_Sequences
#OUTPUT_DIR=output_all_strains_new
#OUTPUT_DIR=output_2ormore_strains_new
#OUTPUT_DIR=output_unique_strains_new

OUTPUT_DIR=core_genes  # OGFilter_max.py
#OUTPUT_DIR=shell_genes # OGFilter_min_max.py
#OUTPUT_DIR=cloud_genes  # OGFilter_min.py

# output_all_strains - core genes - one unique gene per species in all species.
# OK - this threshold works
 MIN_SPECIES=0.92
 MAX_COPIES=50000

# OUTPUT_DIR=output_2ormore_strains
# output_2ormore_strains - shell genes - one unique gene per species in 2 or more species
#MIN_SPECIES=0.15
#MAX_COPIES=50000

# OUTPUT_DIR=output_unique_strains
# output_unique_strains - unique genes - one unique gene per species in one species - required to modify the OGfilter.py script (i.e change greater sign to less than sign)
#MIN_SPECIES=0.15
#MAX_COPIES=50000

python3 OGFilter_max.py -g ${GENECOUNT} -s ${ORTHO_SEQ} -o ${OUTPUT_DIR} -min_species ${MIN_SPECIES} -max_copies ${MAX_COPIES}
