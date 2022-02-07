#!/bin/bash

########################################################################

export assembly=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathRepeatMasker=${path}/results/repeats/${assembly}/RepeatMasker
export pathScripts=${path}/scripts/repeat_annotation

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export prefix=final.contigs
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export prefix=genome_sequence_renamed
fi

#########################################################################

if [ ${assembly} = "NCBI" ]; then
    export prefix=${sp}
fi

########################################################################

if [ ${assembly} = "Ensembl103" ]; then
    export file=`ls ${pathGenomeSequences} | grep ${sp} | grep fa`
    export prefix=`basename ${file} .fa`
fi

#########################################################################

perl ${pathScripts}/combine.repeat.annotations.pl --pathAnnotation1=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation2=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/CombinedRepeatAnnotations.txt --pathOutputBED=${pathRepeatMasker}/CombinedRepeatAnnotations.bed

#########################################################################
