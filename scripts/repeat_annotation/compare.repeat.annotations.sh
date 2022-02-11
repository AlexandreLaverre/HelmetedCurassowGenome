#!/bin/bash

########################################################################

export sp=$1
export method=$2
export assembly=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathRepeatMasker=${path}/results/repeats/${sp}/${assembly}/RepeatMasker
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

perl ${pathScripts}/compare.repeat.annotations.pl --pathAnnotation1=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation2=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/Comparison_Dfam_RepeatModeler.txt

perl ${pathScripts}/compare.repeat.annotations.pl --pathAnnotation2=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation1=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/Comparison_RepeatModeler_Dfam.txt

#########################################################################
