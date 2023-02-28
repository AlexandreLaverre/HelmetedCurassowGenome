#!/bin/bash

export cluster=$1

###########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/gene_family_evolution/all_species/birds/

###########################################################################

cafe5 -i ${pathResults}/gene_families.txt -t ${pathResults}/species_tree.txt -y ${pathResults}/species_tree_annotated.txt -o ${pathResults}/cafe5_results_helmeted_species

###########################################################################
