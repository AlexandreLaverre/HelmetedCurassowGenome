#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
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

fastqc -o ${pathWGS} ${pathFastQ}

#########################################################################
