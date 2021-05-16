#!/bin/bash

export assembly=$1
export annot=$2
export cluster=$3

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/bgeefs/data/necsulea/HelmetedCurassowGenome
else
    if [ ${cluster} = "cloud" ]; then
	export path=/mnt/mydatalocal/HelmetedCurassowGenome
    else
	echo "unknown cluster"
    fi
fi

##############################################################

export pathRepeatMasker=${path}/results/repeats/${assembly}/RepeatMasker/CombinedRepeatAnnotations.txt 
export pathResults=${path}/results/genome_annotation/${assembly}/${annot}
export pathScripts=${path}/scripts/gene_annotation_quality

export release=94

##############################################################

if [ ${annot} = "BRAKER_Ensembl103_multithread" ]; then
    export pathGTF=${pathResults}/braker.gtf  
fi

##############################################################

if [ ${annot} = "GeMoMa" ]; then
    export pathResults=${path}/results/genome_annotation/${assembly}/${annot}/combined
    export pathGTF=${pathResults}/filtered_predictions.gtf
fi

##############################################################

perl ${pathScripts}/overlap.repeats.pl --pathAnnotGTF=${pathGTF} --pathRepeatMasker=${pathRepeatMasker} --pathOutput=${pathResults}/overlap_repeats.txt

##############################################################
