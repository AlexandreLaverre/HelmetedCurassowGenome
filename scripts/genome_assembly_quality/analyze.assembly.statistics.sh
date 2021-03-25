#!/bin/bash

########################################################################

export method=$1	# MEGAHIT or MEGAHIT_RAGOUT
export cluster=$2	# pbil or cloud

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_assembly/${method}
export pathScripts=${path}/scripts/genome_assembly_quality

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathResults}/final.contigs.fa
fi

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathResults}/genome_sequence.fa
fi


#########################################################################

perl ${pathScripts}/analyze.assembly.statistics.pl --pathAssembly=${pathAssembly} --pathOutputStatistics=${pathResults}/assembly.stats.out

#########################################################################
