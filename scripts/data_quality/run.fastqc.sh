#!/bin/bash

export species=$1
export datatype=$2
export cluster=$3
export nthreads=$4

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathData=${path}/data/${datatype}/${species}
export pathScripts=${path}/scripts/data_quality

#########################################################################

export pathFastQ=""

for file in `ls ${pathData} | grep fastq.gz `
do
    export pathFastQ="${pathData}/${file} ${pathFastQ}"
done

#########################################################################

if [ -e ${pathData}/fastqc ]; then
    echo "output dir already there"
else
    mkdir ${pathData}/fastqc
fi

#########################################################################

fastqc -o ${pathData}/fastqc --threads ${nthreads} ${pathFastQ}

#########################################################################
