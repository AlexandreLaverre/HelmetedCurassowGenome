#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/IPLOSS
fi

export pathHAL=${path}/results/genome_assembly/MEGAHIT_RAGOUT

#########################################################################

hal2mafMP.py ${pathHAL}/alignment.hal ${pathHAL}/mafs_by_chr/alignment.maf --numProc 12 --splitBySequence 

#########################################################################
