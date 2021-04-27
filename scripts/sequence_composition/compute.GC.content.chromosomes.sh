#!/bin/bash

########################################################################

export method=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${method}
export pathScripts=${path}/scripts/sequence_composition

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

perl ${pathScripts}/compute.GC.content.chromosomes.pl --pathGenomeSequence=${pathAssembly} --pathOutput=${pathGenomeAssembly}/GCContent_ChromosomeSize.txt

#########################################################################
