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
export pathGenomes=${path}/data/genome_sequences/${source}
export pathAnnotations=${path}/data/genome_annotations/${source}
export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
    export suffix=genome_sequence
fi

#########################################################################

gffread -E ${pathResults}/filtered_predictions.gff -T -o ${pathResults}/filtered_predictions.gtf

gffread -S -y ${pathResults}/filtered_predictions.faa -g ${pathAssembly} ${pathResults}/filtered_predictions.gtf

#########################################################################

perl ${pathScripts}/filter.GeMoMa.predictions.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --minProteinLength=100 --source=GeMoMa --pathOutputGTF=${pathResults}/final_predictions.gtf --pathOutputFasta=${pathResults}/final_predictions.faa

#########################################################################
