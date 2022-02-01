#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export lib=$3
export cluster=$4
export nthreads=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathRepeatModeler=${path}/results/repeats/${sp}/${assembly}/RepeatModeler
export pathResults=${path}/results/repeats/${sp}/${assembly}/RepeatMasker/${lib}
export pathScripts=${path}/scripts/repeat_annotation

## RepeatMasker 4.1.2
## CONS-Dfam 3.3 (Avril 2021)

#########################################################################

if [ -e ${pathResults} ]; then
    echo "path results already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

if [ ${lib} = "Dfam" ]; then
    RepeatMasker -e rmblast -pa ${nthreads} -s -dir ${pathResults} -gff ${pathAssembly}
fi

#########################################################################

if [ ${lib} = "RepeatModeler" ]; then
    RepeatMasker -e rmblast -pa ${nthreads} -s -dir ${pathResults} -gff ${pathAssembly} -lib ${pathRepeatModeler}/repeat_modeler_db-families.fa
fi

#########################################################################
