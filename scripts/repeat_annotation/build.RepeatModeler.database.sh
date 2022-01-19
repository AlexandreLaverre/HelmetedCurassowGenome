#!/bin/bash

set -e 

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

export pathEnsembl=${path}/data/genome_sequences/Ensembl103
export pathNCBI=${path}/data/genome_sequences/NCBI
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathResults=${path}/results/repeats/${sp}/${assembly}/RepeatModeler
export pathScripts=${path}/scripts/repeat_annotation

## RepeatModeler 2.0.1
## RECON 1.08
## RepeatScout 1.0.6
## TRF 4.09
## rmblast 2.11.0

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

if [ ${assembly} = "Ensembl" ]||[ ${assembly} = "Ensembl103" ]; then
    export file=`ls ${pathEnsembl} | grep ${sp}. | grep fa`
    export pathAssembly=${pathEnsembl}/${file}
fi

#########################################################################

if [ ${assembly} = "NCBI" ]; then
    export file=`ls ${pathEnsembl} | grep ${sp}. | grep fa`
    export pathAssembly=${pathNCBI}/${file}
fi

#########################################################################

if [ -e ${pathAssembly} ]; then
    echo "OK, found genome sequence file"
else
    echo "cannot find genome sequence!!"
    exit
fi

#########################################################################

if [ -e ${pathResults} ]; then
    echo "outdir already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

BuildDatabase -name ${pathResults}/repeat_modeler_db ${pathAssembly}

#########################################################################
