#!/bin/bash

export cluster=$1
export minid=$2
export minquerycov=$3
export minsubjectcov=$4

export method="MEGAHIT_RAGOUT"
export brakerset="BRAKER_Ensembl103_multithread"
export genome="genome_sequence_renamed"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathProteins=${path}/data/protein_sequences
export pathResults=${path}/results/genome_annotation/${method}/${brakerset}

##########################################################################

cut -f 1  ${pathResults}/diamond_results/*_minid${minid}_minquerycov${minquerycov}_minsubjectcov${minsubjectcov}.diamond.blastp.out | sort -u > proteins_with_hits_minid${minid}_minquerycov${minquerycov}_minsubjectcov${minsubjectcov}.txt

##########################################################################
