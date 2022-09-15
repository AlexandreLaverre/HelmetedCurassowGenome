#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools
fi

export pathAssembly=${path}/data/genome_sequences/NCBI/${sp}.fa
export pathAnnot=${path}/data/genome_annotation/NCBI/
export pathScripts=${path}/scripts/gene_annotation_quality

#########################################################################

gffread -S -y ${pathAnnot}/${sp}.faa -g ${pathAssembly} ${pathAnnot}/${sp}.gff

gffread -S -x ${pathAnnot}/${sp}.cds.fa -g ${pathAssembly} ${pathAnnot}/${sp}.gff

#########################################################################
