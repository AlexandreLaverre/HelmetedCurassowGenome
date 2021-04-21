#!/bin/bash

########################################################################

export method=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${method}
export pathResults=${path}/results/repeats/${method}/RepeatModeler
export pathScripts=${path}/scripts/repeat_annotation

#########################################################################

RepeatModeler -database ${pathResults}/repeat_modeler_db -pa 8 -LTRStruct >& ${pathResults}/RepeatModeler.out 

#########################################################################
