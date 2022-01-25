#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3
export nthreads=$4

########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/RNASeq_alignments/${sp}/${assembly}
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}

########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export prefix="genome_sequence_renamed"
    export pathGenomeSequence=${pathGenomeAssembly}/${assembly}/${prefix}.fa
fi

########################################################################

if [ -e ${pathResults} ]; then
    echo "path results exists"
else
    mkdir -p ${pathResults}
fi

########################################################################

hisat2-build --seed 19 -p ${nthreads} ${pathGenomeSequence} ${pathResults}/${prefix}

########################################################################
