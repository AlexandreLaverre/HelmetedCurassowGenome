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
export pathAnnotScripts=${path}/scripts/gene_annotation

#########################################################################

perl ${pathScripts}/rename.GeMoMa.transcripts.pl --pathInputGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25_maxXStop0.1.gtf --pathOutputGTF=${pathResults}/filtered_GeMoMa_annotations.gtf

cp ${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25_maxXStop0.1.cds.fa ${pathResults}/filtered_GeMoMa_annotations.cds.fa

cp ${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25_maxXStop0.1.faa ${pathResults}/filtered_GeMoMa_annotations.faa

#########################################################################

## add gene id and transcript id in fasta file

perl ${pathAnnotScripts}/format.GeMoMa.sequences.pl --pathAnnotGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25_maxXStop0.1.gtf  --pathSequences=${pathResults}/filtered_GeMoMa_annotations.faa --pathOutput=${pathResults}/filtered_GeMoMa_annotations_formatted.faa

perl ${pathAnnotScripts}/format.GeMoMa.sequences.pl --pathAnnotGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25_maxXStop0.1.gtf  --pathSequences=${pathResults}/filtered_GeMoMa_annotations.cds.fa --pathOutput=${pathResults}/filtered_GeMoMa_annotations_formatted.cds.fa

#############################################################################
