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
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeSequences=${path}/data/genome_sequences/${assembly}
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

if [ -e ${pathRepeatMasker}/Comparison_Dfam_RepeatModeler.txt ]; then
    echo "Dfam vs RepeatModeler already done"
else
    perl ${pathScripts}/compare.repeat.annotations.pl --pathAnnotation1=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation2=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/Comparison_Dfam_RepeatModeler.txt
fi

#########################################################################

if [ -e ${pathRepeatMasker}/Comparison_RepeatModeler_Dfam.txt ]; then
    echo "RepeatModeler vs Dfam already done"
else
    perl ${pathScripts}/compare.repeat.annotations.pl --pathAnnotation2=${pathRepeatMasker}/Dfam/${prefix}.fa.out --pathAnnotation1=${pathRepeatMasker}/RepeatModeler/${prefix}.fa.out --pathOutput=${pathRepeatMasker}/Comparison_RepeatModeler_Dfam.txt
fi
#########################################################################
