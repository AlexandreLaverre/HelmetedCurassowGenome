#!/bin/bash

########################################################################

export sp=$1
export method=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${method}
export pathScripts=${path}/scripts/genome_assembly_quality

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence.fa
    export suffix=genome_sequence
fi

#########################################################################

perl ${pathScripts}/rename.chromosomes.pl --pathGenomeSequence=${pathAssembly} --pathChromosomeCorrespondence=${pathGenomeAssembly}/sequence_names.txt --pathOutput=${pathGenomeAssembly}/genome_sequence_renamed.fa 

#########################################################################
