#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeSequences=${path}/data/genome_sequences/${assembly}
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
    export file=`ls ${pathGenomeSequences} | grep ${sp} | grep fa.gz | grep -v sm`

    if [ -e ${pathGenomeSequences}/${file} ]; then
	export prefix=`basename ${file} .fa.gz`
    else
	export file=`ls ${pathGenomeSequences} | grep ${sp} | grep fa | grep -v sm`

	if [ -e ${pathGenomeSequences}/${file} ]; then
	    export prefix=`basename ${file} .fa`
	else
	    echo "cannot find genome sequence file"
	    exit
	fi
    fi
fi

#########################################################################

if [ -e ${pathRepeatMasker}/CombinedRepeatAnnotations.txt ]&& [ -e ${pathRepeatMasker}/CombinedRepeatAnnotations.bed ] ; then
    echo "already done"
else
    perl ${pathScripts}/combine.repeat.annotations.pl --pathAnnotation1=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation2=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/CombinedRepeatAnnotations.txt --pathOutputBED=${pathRepeatMasker}/CombinedRepeatAnnotations.bed
fi
#########################################################################
