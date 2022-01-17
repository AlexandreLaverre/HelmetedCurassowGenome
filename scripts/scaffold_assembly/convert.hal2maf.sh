#!/bin/bash

export sp=$1
export cluster=$2
export nthreads=$3

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${path}/results/genome_assembly/${sp}/MEGAHIT_RAGOUT

#########################################################################

hal2mafMP.py ${pathHAL}/alignment.hal ${pathHAL}/mafs_by_chr/alignment.maf --numProc ${nthreads} --splitBySequence

#########################################################################
