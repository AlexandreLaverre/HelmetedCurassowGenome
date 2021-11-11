#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathIndex=${path}/data/genome_indexes/MEGAHIT_RAGOUT

#########################################################################

bwa index ${pathIndex}/genome_sequence_renamed.fa

#########################################################################
