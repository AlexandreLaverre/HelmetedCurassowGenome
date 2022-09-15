#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
    export pathTools=/beegfs/home/${USER}/Tools/OrthoFinder/tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathCDS=${path}/data/coding_sequences/Ensembl103
export pathScripts=${path}/scripts/gene_annotation_quality

#########################################################################

export cdsfile=`ls ${pathCDS} | grep "${sp}\." | grep fa`

#########################################################################

python ${pathTools}/primary_transcript.py ${pathCDS}/${cdsfile}

#########################################################################

