#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_assembly/MEGAHIT_RAGOUT

#########################################################################

ragout -o ${pathResults} -s maf --refine -t 24 ragout.recipe.${cluster} 

#########################################################################
