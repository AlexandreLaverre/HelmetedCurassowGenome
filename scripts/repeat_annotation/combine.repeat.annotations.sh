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
export pathRepeatMasker=${path}/results/repeats/${method}/RepeatMasker
export pathScripts=${path}/scripts/repeat_annotation

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export prefix=final.contigs
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export prefix=genome_sequence_renamed
fi

#########################################################################

perl ${pathScripts}/combine.repeat.annotations.pl --pathAnnotation1=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation2=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/CombinedRepeatAnnotations.txt

cut -f 1-3 ${pathRepeatMasker}/CombinedRepeatAnnotations.txt | sed '1d' > ${pathRepeatMasker}/CombinedRepeatAnnotations.bed

#########################################################################

