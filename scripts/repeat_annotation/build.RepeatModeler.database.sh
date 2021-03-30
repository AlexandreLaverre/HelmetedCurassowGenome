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
export pathResults=${path}/results/repeats/${method}
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

BuildDatabase -name ${pathResults}/repeat_modeler_db ${pathAssembly}

#########################################################################
