#!/bin/bash

export species=$1
export cluster=$2
export nthreads=$3

#############################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathRNASeq=${path}/data/RNASeq/${species}
export pathDocs=${path}/docs
export pathResults=${path}/results/genome_assembly/${species}/MEGAHIT
export pathScripts=${path}/scripts/scaffold_assembly

export pathIndex=${pathResults}/final.contigs

#############################################################################

samtools merge -@ ${nthreads} -o ${pathResults}/accepted_hits_all_samples.bam ${pathResults}/hisat*/accepted_hits.bam

#############################################################################
