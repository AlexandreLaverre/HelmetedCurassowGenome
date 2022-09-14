#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/filter_gene_annotation

#########################################################################

export pathDiamondResults=${pathResults}/diamond_results/SignificantHits_MinProteinFraction0.25_MaxEvalue0.001_MaxGapGraction0.1.txt

#########################################################################

## quantile(chicken protein lengths, p=0.01) = 68 

perl ${pathScripts}/filter.GeMoMa.predictions.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --pathDiamondResults=${pathDiamondResults} --minProteinLength=70 --pathOverlapRepeats=${pathResults}/overlap_repeats.txt --maxFractionRepeats=0.25 --source=GeMoMa --pathOutputGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25.gtf --pathOutputFasta=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25.faa 

#########################################################################
