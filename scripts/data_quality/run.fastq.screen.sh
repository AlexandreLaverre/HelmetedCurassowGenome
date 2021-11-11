#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathWGS=${path}/data/WGS
export pathScripts=${path}/scripts/data_quality

#########################################################################

export pathFastQ=""

for file in `ls ${pathWGS} | grep fastq.gz `
do
    export pathFastQ="${pathWGS}/${file} ${pathFastQ}"
done

#########################################################################

if [ -e ${pathWGS}/fastq_screen ]; then
    echo "output dir already there"
else
    mkdir ${pathWGS}/fastq_screen
fi

#########################################################################

fastq_screen --aligner bowtie2  --subset 1000000 --outdir ${pathWGS}/fastq_screen ${pathFastQ}

#########################################################################
