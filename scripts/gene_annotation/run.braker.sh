#!/bin/bash

########################################################################

export method=$1
export lib=$2
export cluster=$3
export cores=$4

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${method}
export pathResults=${path}/results/genome_annotation/${method}/BRAKER_${lib}
export pathProteins=${path}/data/protein_sequences/${lib}
export pathScripts=${path}/scripts/genome_assembly_quality

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
    export suffix=genome_sequence
fi

#########################################################################

braker.pl --genome=${pathAssembly} --prot_seq=${pathProteins}/AllProteins.fa --cores=${cores} --epmode --softmasking

if [ -e ${pathResults} ]; then
    echo "output dir already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################
