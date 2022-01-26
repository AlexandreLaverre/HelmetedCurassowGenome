#!/bin/bash

export cluster=$1
export annot=$2
export type=$3
export threads=$4

export assembly="MEGAHIT_RAGOUT"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools/OrthoFinder/tools
fi

export pathResults=${path}/results/genome_annotation/${assembly}/${annot}

##########################################################################

ulimit -n 50000

## OrthoFinder 2.5.4
## MAFFT v7.453
## fasttree ubuntu package 2.1.11-1

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
