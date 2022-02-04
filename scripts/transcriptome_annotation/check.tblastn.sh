#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export refsp=$3
export source=$4
export cluster=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathProteinSequences=${path}/data/protein_sequences/${source}/primary_transcripts
export pathTranscriptomeAssembly=${path}/results/transcriptome_assembly/${target}
export pathResults=${path}/results/transcriptome_assembly/${target}/tblastn_results
export pathScripts=${path}/scripts/transcriptome_annotation

########################################################################

export protfile=`ls ${pathProteinSequences} | grep ${refsp} | grep fa`

if [ -e ${pathProteinSequences}/${protfile} ]; then
    echo "using as an input "${pathProteinSequences}/${protfile}
else
    echo "cannot find protein sequence file for "${refsp}
    exit
fi

#########################################################################

if [ -e ${pathResults}/${refsp}_vs_${suffix}.tblastn.out ]; then

    export last=`tail -n 1 ${pathResults}/${refsp}_vs_${suffix}.tblastn.out | cut -f 1`
    export tot=`grep -c ">" ${pathProteinSequences}/${protfile}`
    export index=`grep ">" ${pathProteinSequences}/${protfile} | grep -n ${last} | cut -f 1 -d ':'`

    export ratio=$(($index * 100 / $tot))

    echo "index "${index}" out of "${tot}" "${ratio}"% done";
else
    echo "tblastn not yet done"
fi

#########################################################################
