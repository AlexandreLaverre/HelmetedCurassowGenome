#!/bin/bash

export type=$1  ## tree inference method
export cluster=$2

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

export pathResults=${path}/results/gene_families/OrthoFinder/all_species/${type}
export pathScripts=${path}/scripts/gene_families

##########################################################################

export dir=`ls ${pathResults} | grep Results`
export pathOrthoGroups=${pathResults}/${dir}/Phylogenetic_Hierarchical_Orthogroups

##########################################################################

perl ${pathScripts}/extract.orthogroup.correspondence.pl --pathN0=${pathOrthoGroups}/N0.tsv --pathN1=${pathOrthoGroups}/N1.tsv --pathN2=${pathOrthoGroups}/N2.tsv --pathOutput=${pathOrthoGroups}/N0N1N2_correspondence.txt

##########################################################################
