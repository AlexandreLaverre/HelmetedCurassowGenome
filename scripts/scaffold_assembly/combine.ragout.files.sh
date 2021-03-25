#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_assembly/MEGAHIT_RAGOUT

#########################################################################

cat ${pathResults}/Pauxi_pauxi_scaffolds.fasta ${pathResults}/Pauxi_pauxi_unplaced.fasta > ${pathResults}/genome_sequence.fa
 
#########################################################################
