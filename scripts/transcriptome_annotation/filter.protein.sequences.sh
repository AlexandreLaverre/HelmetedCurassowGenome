#!/bin/bash

export cluster=$1
export annot=$2

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools/OrthoFinder/tools
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome/
    export pathTools=/beegfs/home/necsulea/Tools/OrthoFinder/tools
fi

export pathSequences=${path}/data/protein_sequences/${annot}

##########################################################################

for file in `ls ${pathSequences} | grep fa`
do
    if [ -e ${pathSequences}/primary_transcripts/${file} ]; then
	echo "already extracted primary transcripts"
    else
	echo "extracting primary transcripts for "${file}
	python ${pathTools}/primary_transcript.py ${pathSequences}/${file}
    fi
done

##########################################################################
