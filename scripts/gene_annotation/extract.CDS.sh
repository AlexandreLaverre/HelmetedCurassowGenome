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

gffread -S -x ${pathResults}/final_annotations.cds.fa -g ${pathAssembly} ${pathResults}/final_annotations.gtf

perl ${pathScripts}/format.GeMoMa.proteins.pl --pathAnnotGTF=${pathResults}/final_annotations.gtf --pathProteins=${pathResults}/final_annotations.cds.fa --pathOutput=${pathResults}/final_annotations_formatted.cds.fa

#########################################################################
