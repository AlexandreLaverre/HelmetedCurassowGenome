#!/bin/bash

export cluster=$1

export assembly="MEGAHIT_RAGOUT"
export brakerset="Ensembl103_multithread"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_annotation/${assembly}/BRAKER_${brakerset}

##########################################################################

grep intron ${pathResults}/braker.gtf  | cut -f 9 | cut -f 4 -d ' ' | cut -f 2 -d '"' | sort -u > ${pathResults}/multiexonic_genes.txt

grep intron ${pathResults}/braker.gtf  | cut -f 9 | cut -f 2 -d ' ' | cut -f 2 -d '"' | sort -u > ${pathResults}/multiexonic_genes_proteins.txt

##########################################################################

