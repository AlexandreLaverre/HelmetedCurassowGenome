#!/bin/bash

set -e

export cluster=$1

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

export pathChainFiles=${path}/data/liftOver_files
export pathResults=${path}/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse
export pathScripts=${path}/scripts/noncoding_element_evolution

##########################################################################

liftOver -minMatch=0.25 -multiple ${pathResults}/combined_peaks_mm10.bed ${pathChainFiles}/mm10ToGalGal6.over.chain.gz ${pathResults}/combined_peaks_galGal6_multiple.bed ${pathResults}/combined_peaks_galGal6_unmapped.bed 

##########################################################################

liftOverMerge -mergeGap=100 ${pathResults}/combined_peaks_galGal6_multiple.bed ${pathResults}/combined_peaks_galGal6_merged100bp.bed

##########################################################################
