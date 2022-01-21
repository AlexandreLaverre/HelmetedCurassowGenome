#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathScripts=${path}/scripts/sequence_composition

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

perl ${pathScripts}/compute.GC.content.chromosomes.pl --pathGenomeSequence=${pathAssembly} --pathOutput=${pathGenomeAssembly}/GCContent_ChromosomeSize.txt

#########################################################################
