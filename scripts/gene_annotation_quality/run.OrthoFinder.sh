#!/bin/bash

export cluster=$1
export annot=$2
export type=$3

export assembly="MEGAHIT_RAGOUT"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathResults=${path}/results/genome_annotation/${assembly}/${annot}

##########################################################################

ulimit -n 50000

if [ ${type} = "final" ]; then
    orthofinder -f ${pathResults}/OrthoFinder_${type} -t 30 -I 2 -M msa 
fi

##########################################################################
