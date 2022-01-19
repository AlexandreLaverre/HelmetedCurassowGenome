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
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${method}
export pathResults=${path}/results/repeats/${method}
export pathScripts=${path}/scripts/repeat_annotation

## RepeatModeler 2.0.1
## RECON 1.08
## RepeatScout 1.0.6
## TRF 4.09
## rmblast 2.11.0

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

BuildDatabase -name ${pathResults}/repeat_modeler_db ${pathAssembly}

#########################################################################
