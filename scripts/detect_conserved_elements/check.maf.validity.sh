#!/bin/bash

export dataset=$1
export refsp=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

export pathAln=${path}/results/whole_genome_alignments/${dataset}/${refsp}
export pathScripts=${path}/scripts/detect_conserved_elements

#########################################################################

for file in `ls ${pathAln} | grep maf`
do
    echo ${file} 
    perl ${pathScripts}/check.maf.validity.pl --pathMAF=${pathAln}/${file} 
done

#########################################################################
