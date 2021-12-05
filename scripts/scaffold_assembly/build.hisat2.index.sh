#!/bin/bash

export sp=$1
export cluster=$2
export threads=$3

##############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathMEGAHIT=${path}/results/genome_assembly/${sp}/MEGAHIT

##############################################################################

hisat2-build -p ${threads} ${pathMEGAHIT}/final.contigs.fa ${pathMEGAHIT}/final.contigs

##############################################################################
