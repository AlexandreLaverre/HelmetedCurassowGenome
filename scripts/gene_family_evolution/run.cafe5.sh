#!/bin/bash

export cluster=$1

###########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome/all_species/birds/
fi

export pathResults=${path}/results/gene_family_evolution

###########################################################################

cafe5 -i ${pathResults}/gene_families.txt -t ${pathResults}/species_tree_annotated.txt -y ${pathResults}/cafe5_results_helmeted_species.txt

###########################################################################
