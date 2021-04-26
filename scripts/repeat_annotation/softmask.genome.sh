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
export pathResults=${path}/results/repeats/${method}/RepeatMasker
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
