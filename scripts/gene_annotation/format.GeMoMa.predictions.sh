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
    export pathTools=/mnt/mydatalocal/Tools
fi

export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
fi

#########################################################################

if [ ${assembly} = "NCBI" ]; then
    export pathGenomeSequences=${path}/data/genome_sequences/${assembly}
    export pathAssembly=${pathGenomeSequences}/${sp}.fa
fi

#########################################################################

gffread -E ${pathResults}/filtered_predictions.gff -T -o ${pathResults}/filtered_predictions.gtf

gffread -S -y ${pathResults}/filtered_predictions.faa -g ${pathAssembly} ${pathResults}/filtered_predictions.gff

gffread -S -x ${pathResults}/filtered_predictions.cds.fa -g ${pathAssembly} ${pathResults}/filtered_predictions.gff

#########################################################################

perl ${pathScripts}/format.GeMoMa.proteins.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --pathOutput=${pathResults}/filtered_predictions_formatted.faa

#########################################################################
