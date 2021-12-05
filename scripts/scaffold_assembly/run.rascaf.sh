#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=//home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathMEGAHIT=${path}/results/genome_assembly/${sp}/MEGAHIT
export pathResults=${path}/results/genome_assembly/${sp}/rascaf

#########################################################################

if [ -e ${pathResults} ]; then
    echo "dir results already there"
else
    mkdir ${pathResults}
fi

#########################################################################

rascaf -o ${pathResults}/final -b ${pathMEGAHIT}/accepted_hits_all_samples.bam -f ${pathMEGAHIT}/final.contigs.fa -ms 5 

#########################################################################
