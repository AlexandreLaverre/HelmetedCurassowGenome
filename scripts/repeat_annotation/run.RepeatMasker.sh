#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export lib=$3
export cluster=$4
export nthreads=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeSequences=${path}/data/genome_sequences/${assembly}
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathRepeatModeler=${path}/results/repeats/${sp}/${assembly}/RepeatModeler
export pathResults=${path}/results/repeats/${sp}/${assembly}/RepeatMasker/${lib}
export pathScripts=${path}/scripts/repeat_annotation

## RepeatMasker 4.1.2-p1
## CONS-Dfam 3.3 (Avril 2021)

##Using Master RepeatMasker Database: /ifb/data/mydatalocal/Tools/RepeatMasker/Libraries/RepeatMaskerLib.h5
## Title    : Dfam
## Version  : 3.3
## Date     : 2020-11-09
##Families : 273,693


# cd-hit-v4.8.1-2019-0228
# LTR_retriever-2.9.0
# mafft-7.475-without-extensions -> mafft 7.453
# RECON-1.08
# RepeatModeler-2.0.1
# rmblast-2.11.0
# trf409.linux64
# genometools-1.6.1
# NINJA-0.95-cluster_only
# RepeatMasker
# RepeatScout-1.0.6
# salmon-1.6.0_linux_x86_64
# trinityrnaseq-v2.13.2
# jellyfish 2.3.0

#########################################################################

if [ -e ${pathResults} ]; then
    echo "path results already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
fi

#########################################################################

if [ ${assembly} = "NCBI" ]; then
    export pathAssembly=${pathGenomeSequences}/${sp}.fa
fi

#########################################################################

if [ ${assembly} = "Ensembl103" ]; then
    export file=`ls ${pathGenomeSequences} | grep ${sp} | grep fa`
    export pathAssembly=${pathGenomeSequences}/${file}
fi

#########################################################################

if [ ${lib} = "Dfam" ]; then
    RepeatMasker -e rmblast -pa ${nthreads} -s -dir ${pathResults} -gff ${pathAssembly}
fi

#########################################################################

if [ ${lib} = "RepeatModeler" ]; then
    RepeatMasker -e rmblast -pa ${nthreads} -s -dir ${pathResults} -gff ${pathAssembly} -lib ${pathRepeatModeler}/repeat_modeler_db-families.fa
fi

#########################################################################
