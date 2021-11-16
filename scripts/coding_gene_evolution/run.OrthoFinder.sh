#!/bin/bash

export type=$1  ## tree inference method
export cluster=$2
export threads=$3

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathResults=${path}/results/coding_gene_evolution/

##########################################################################

ulimit -n 50000

## OrthoFinder 2.5.4
## MAFFT v7.453
## fasttree ubuntu package 2.1.11-1
## IQ-TREE multicore version 1.6.12 for Linux 64-bit built Mar 23 2020

##########################################################################

if [ ${type} = "fasttree" ]; then
    orthofinder -f ${pathResults} -o ${pathResults}/OrthoFinder_${type} -t ${threads} -I 2 -M msa -y -T ${type}
fi

##########################################################################

## we use previously inferred species tree

if [ ${type} = "iqtree" ]; then
    orthofinder -f ${pathResults} -o ${pathResults}/OrthoFinder_${type} -t ${threads} -I 2 -M msa -y -T ${type} -s ${pathResults}/species_tree.txt
fi

##########################################################################
