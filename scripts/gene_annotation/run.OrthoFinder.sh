#!/bin/bash

export cluster=$1
export annot=$2
export type=$3
export threads=$4

export assembly="MEGAHIT_RAGOUT"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathResults=${path}/results/genome_annotation/${assembly}/${annot}

##########################################################################

ulimit -n 50000

##########################################################################

## to determine homologous gene families for annotation

if [ ${type} = "final" ]; then
    orthofinder -f ${pathResults}/OrthoFinder_${type} -t ${threads} -I 2 -M msa
fi

##########################################################################

## to roughly cluster protein sequences

if [ ${type} = "filtered" ]; then
    orthofinder -f ${pathResults}/OrthoFinder_${type} -t ${threads} -I 2
fi

##########################################################################
