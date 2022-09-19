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

export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/filter_gene_annotation

#########################################################################

perl ${pathScripts}/rename.GeMoMa.transcripts.pl --pathInputGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25_maxXStop0.1.gtf --pathOutputGTF=${pathResults}/filtered_GeMoMa_annotations.gtf

#########################################################################
