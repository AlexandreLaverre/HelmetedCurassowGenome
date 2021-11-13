#!/bin/bash

export assembly=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/combined
export pathOrthoFinder=${pathResults}/OrthoFinder_filtered/OrthoFinder
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

export pathOrthoGroups=`ls ${pathOrthoFinder}/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv `

#########################################################################

perl ${pathScripts}/filter.GeMoMa.predictions.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --pathOrthoGroups=${pathOrthoGroups} --minProteinLength=100 --pathOverlapRepeats=${pathResults}/overlap_repeats.txt --maxFractionRepeats=0.5 --source=GeMoMa --pathOutputGTF=${pathResults}/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.gtf --pathOutputFasta=${pathResults}/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa 

#########################################################################
