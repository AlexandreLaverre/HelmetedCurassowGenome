#!/bin/bash

export sp=$1
export cluster=$2
export nthreads=$3

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${path}/results/genome_assembly/${sp}/MEGAHIT_RAGOUT

#########################################################################

if [ -e ${pathHAL}/mafs_by_chr ]; then
    echo "dir output exists"
else
    mkdir ${pathHAL}/mafs_by_chr
fi

#########################################################################

docker run -v ${path}:/ifb/data/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 hal2mafMP.py ${pathHAL}/alignment.hal ${pathHAL}/mafs_by_chr/alignment.maf --numProc ${nthreads} --splitBySequence

#########################################################################

