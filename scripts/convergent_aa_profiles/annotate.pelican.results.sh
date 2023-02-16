#!/bin/bash

set -e

export spset=$1
export dataset=$2
export run=$3
export cluster=$4

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

export pathResults=${path}/results/coding_gene_evolution/${spset}/${dataset}/pelican_output_${run}
export pathOrthoFinder=${path}/results/gene_families/OrthoFinder/all_species/iqtree
export dir=`ls ${pathOrthoFinder} | grep Results`
export pathOrthoGroups=${pathOrthoFinder}/${dir}/Phylogenetic_Hierarchical_Orthogroups
export pathScripts=${path}/scripts/convergent_aa_profiles

##########################################################################

if [ ${dataset} = "birds" ]; then
    export orthogroup="N2"
fi

if [ ${dataset} = "squamates" ]; then
    export orthogroup="N1"
fi

##########################################################################

for res in best_sites all_sites
do
    echo "perl ${pathScripts}/annotate.pelican.results.pl --pathPelicanResults=${pathResults}/${res}.tsv --orthogroup=${orthogroup} --pathOrthoGroupAnnotation=${pathOrthoGroups}/${orthogroup}_annotations.tsv --pathOrthoGroupCorrespondence=${pathOrthoGroups}/N0N1N2_correspondence.txt --pathOutput=${pathResults}/${res}_annotated.tsv"
done

##########################################################################
