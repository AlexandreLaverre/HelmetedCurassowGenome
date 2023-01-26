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

export pathDocs=${path}/docs
export pathData=${path}/data/ENCODE_ATAC-seq/Mouse
export pathResults=${path}/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse
export pathScripts=${path}/scripts/noncoding_element_evolution

##########################################################################

if [ -e ${pathResults} ]; then
    echo "output dir already there"
else
    mkdir -p ${pathResults}
fi

##########################################################################

export samples=""
export paths=""

for file in `ls ${pathData} | grep GSM`
do
    export sample=`echo ${file} | cut -f 1 -d '_'`
    export samples=${sample},${samples}
    export paths=${pathData}/${file},${paths}
done

##########################################################################

perl ${pathScripts}/combine.peaks.pl --pathSampleInfo=${pathDocs}/SraRunTable_ENCODE_ATAC_seq_Mus_musculus.tsv --samples=${samples} --pathsCoordinates=${paths} --pathOutput=${pathResults}/combined_peaks_mm10.txt --pathOutputBED=${pathResults}/combined_peaks_mm10.bed

##########################################################################
