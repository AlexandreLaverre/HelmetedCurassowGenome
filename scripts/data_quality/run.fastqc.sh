#!/bin/bash

export sp="HelmetedCurassow"

#########################################################################

export path=/beegfs/data/necsulea/IPLOSS
export pathWGS=${path}/data/WGS/${sp}
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
