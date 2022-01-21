#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_assembly/${sp}/MEGAHIT_RAGOUT

#########################################################################

cat ${pathResults}/${sp}_scaffolds.fasta ${pathResults}/${sp}_unplaced.fasta > ${pathResults}/genome_sequence.fa
 
#########################################################################
