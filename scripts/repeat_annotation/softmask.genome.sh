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
export pathResults=${path}/results/repeats/${sp}/${method}/RepeatMasker
export pathScripts=${path}/scripts/repeat_annotation

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export prefixAssembly=final.contigs
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export prefixAssembly=genome_sequence_renamed
fi

#########################################################################

bedtools maskfasta -soft -fi ${pathGenomeAssembly}/${prefixAssembly}.fa -bed ${pathResults}/CombinedRepeatAnnotations.bed  -fo ${pathGenomeAssembly}/${prefixAssembly}_sm.fa

#########################################################################
