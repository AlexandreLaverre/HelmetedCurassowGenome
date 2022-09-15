#!/bin/bash

export sp=$1
export cluster=$2

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

export pathNCBI=${path}/data/genome_annotations/NCBI
export pathGeMoMa=${path}/results/genome_annotation/${sp}/NCBI/GeMoMa/combined
export pathScripts=${path}/scripts/filter_gene_annotation

#############################################################################

if [ -e ${pathNCBI}/${sp}.gtf ]; then
    echo "gtf file already there"
else
    gffread -E ${pathNCBI}/${sp}.gff -T -o ${pathNCBI}/${sp}.gtf
fi

#############################################################################

for minlen in 70 100
do
    for maxrep in 0.25 0.5
    do

	perl ${pathScripts}/compare.annotations.pl --pathAnnot1=${NCBI}/${sp}.gtf --pathAnnot2=${pathGeMoMa}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.gtf --pathOutput=${pathGeMoMa}/NCBI_vs_filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.txt

	perl ${pathScripts}/compare.annotations.pl --pathAnnot2=${NCBI}/${sp}.gtf --pathAnnot1=${pathGeMoMa}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.gtf --pathOutput=${pathGeMoMa}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}_vs_NCBI.txt
	
    done
done
#############################################################################
