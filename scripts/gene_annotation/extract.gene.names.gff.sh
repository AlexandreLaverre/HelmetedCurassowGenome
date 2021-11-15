#!/bin/bash

export sp=$1
export source=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathAnnot=${path}/data/genome_annotations/${source}
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

export annotfile=`ls ${pathAnnotations} | grep ${ref}'\.' | grep gff`

#########################################################################

perl ${pathScripts}/extract.gene.names.gff.pl --pathGFF=${pathAnnot}/${annotfile} --pathOutput=${pathAnnot}/${ref}.gene.names.txt

#########################################################################
