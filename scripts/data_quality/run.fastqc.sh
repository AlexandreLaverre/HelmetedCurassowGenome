#!/bin/bash

export species=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=//home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathWGS=${path}/data/WGS/${species}
export pathScripts=${path}/scripts/data_quality

#########################################################################

export pathFastQ=""

for file in `ls ${pathWGS} | grep fastq.gz `
do
    export pathFastQ="${pathWGS}/${file} ${pathFastQ}"
done

#########################################################################

if [ -e ${pathWGS}/fastqc ]; then
    echo "output dir already there"
else
    mkdir ${pathWGS}/fastqc
fi

#########################################################################

fastqc -o ${pathWGS}/fastqc ${pathFastQ}

#########################################################################
