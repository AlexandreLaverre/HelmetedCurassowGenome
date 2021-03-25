#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomes=${path}/data/genome_sequences
export pathScripts=${path}/scripts/ragout_assembly

#########################################################################

if [ -e ${pathGenomes}/${sp}/genome_sm.fa ]; then
    perl ${pathScripts}/cleanup.fasta.names.pl --pathInput=${pathGenomes}/${sp}/genome_sm.fa --pathOutput=${pathGenomes}/${sp}/genome_sm_clean.fa
fi

#########################################################################

