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

export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

perl ${pathScripts}/rename.GeMoMa.transcripts.pl --pathInputGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.5.gtf --pathOutputGTF=${pathResults}/final_annotations.gtf

#########################################################################
