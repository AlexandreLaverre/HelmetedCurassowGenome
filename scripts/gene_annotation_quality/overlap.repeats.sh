#!/bin/bash

export sp=$1
export assembly=$2
export annot=$3
export cluster=$4

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
else
    if [ ${cluster} = "cloud" ]; then
	export path=/mnt/mydatalocal/HelmetedCurassowGenome
    else
	echo "unknown cluster"
    fi
fi

##############################################################

export pathRepeatMasker=${path}/results/repeats/${sp}/${assembly}/RepeatMasker/CombinedRepeatAnnotations.txt
export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/${annot}
export pathScripts=${path}/scripts/gene_annotation_quality

##############################################################

if [ ${annot} = "BRAKER_Ensembl103_multithread" ]; then
    export pathGTF=${pathResults}/braker.gtf
fi

##############################################################

if [ ${annot} = "GeMoMa" ]; then
    export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/${annot}/combined
    export pathGTF=${pathResults}/filtered_predictions.gtf
fi

##############################################################

perl ${pathScripts}/overlap.repeats.pl --pathAnnotGTF=${pathGTF} --pathRepeatMasker=${pathRepeatMasker} --pathOutput=${pathResults}/overlap_repeats.txt

##############################################################
