#!/bin/bash

export cluster=$1
export method="MEGAHIT_RAGOUT"
export brakerset="BRAKER_Ensembl103_multithread"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathResults=${path}/results/genome_annotation/${method}/${brakerset}

##########################################################################

orthofinder -f ${pathResults}/OrthoFinder -t 30 

##########################################################################
