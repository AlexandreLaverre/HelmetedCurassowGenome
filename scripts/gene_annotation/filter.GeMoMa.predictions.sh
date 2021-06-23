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

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
fi

#########################################################################

#gffread -E ${pathResults}/filtered_predictions.gff -T -o ${pathResults}/filtered_predictions.gtf

#gffread -S -y ${pathResults}/filtered_predictions.faa -g ${pathAssembly} ${pathResults}/filtered_predictions.gff

#########################################################################

perl ${pathScripts}/filter.GeMoMa.predictions.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --minProteinLength=50 --pathOverlapRepeats=${pathResults}/overlap_repeats.txt --maxFractionRepeats=0.5 --source=GeMoMa --pathOutputGTF=${pathResults}/final_predictions.gtf --pathOutputFasta=${pathResults}/final_predictions.faa --pathOutputFullFasta=${pathResults}/filtered_predictions_formatted.faa

#########################################################################
