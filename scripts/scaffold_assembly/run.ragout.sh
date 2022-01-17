#!/bin/bash

export sp=$1
export cluster=$2
export ncores=$3

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_assembly/${sp}/MEGAHIT_RAGOUT

#########################################################################

ragout -o ${pathResults} -s maf --refine -t ${ncores} ragout.recipe.${cluster}.${sp}

#########################################################################
