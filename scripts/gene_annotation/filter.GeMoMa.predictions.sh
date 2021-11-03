#!/bin/bash

export assembly=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

perl ${pathScripts}/filter.GeMoMa.predictions.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --pathOrthoGroups=NA --minProteinLength=100 --pathOverlapRepeats=${pathResults}/overlap_repeats.txt --maxFractionRepeats=0.5 --source=GeMoMa --pathOutputGTF=${pathResults}/filtered_predictions_minLength100_maxFractionRepeats0.5.gtf --pathOutputFasta=${pathResults}/filtered_predictions_minLength100_maxFractionRepeats0.5.faa --pathOutputFullFasta=${pathResults}/filtered_predictions_formatted.faa

#########################################################################
