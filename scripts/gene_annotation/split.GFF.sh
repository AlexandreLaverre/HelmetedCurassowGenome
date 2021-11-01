#!/bin/bash

export ref=$1
export source=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathAnnotations=${path}/data/genome_annotations/${source}
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

if [ -e ${pathAnnotations}/parts ]; then
    echo "dir output exists"
else
    mkdir -p ${pathAnnotations}/parts
fi

#########################################################################

export annotfile=`ls ${pathAnnotations} | grep ${ref}'\.' | grep gff`
export prefix=`echo ${annotfile} | cut -f 1 -d '-'`

echo "annotation file ".${annotfile}." prefix "${prefix}

#########################################################################

perl ${pathScripts}/split.GFF.pl --pathGFF=${pathAnnotations}/${annotfile} --maxChrSize=50000000 --dirOutput=${pathAnnotations}/parts --prefixOutput=${prefix}

#########################################################################
