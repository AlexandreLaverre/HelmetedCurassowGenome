#!/bin/bash

########################################################################

export method=$1
export lib=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${method}
export pathResults=${path}/results/repeats/${method}/RepeatMasker/${lib}
export pathScripts=${path}/scripts/repeat_annotation

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

if [ ${lib} = "Dfam" ]; then
    RepeatMasker -e rmblast -pa 8 -s -dir ${pathResults} -gff ${pathAssembly}
fi

#########################################################################
