#!/bin/bash

export assembly=$1
export ref=$2
export source=$3
export cluster=$4
export threads=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathGenomes=${path}/data/genome_sequences/${source}
export pathAnnotations=${path}/data/genome_annotations/${source}
export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/${ref}

#########################################################################

if [ -e ${pathResults} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
    export suffix=genome_sequence
fi

#########################################################################

export genomefile=`ls ${pathGenomes} | grep ${ref} | grep fa`
export annotfile=`ls ${pathAnnotations} | grep ${ref} | grep gff`

#########################################################################

if [ -e ${pathResults}/final_annotation.gff ]; then
    eho "already done"
else
    java -jar /mnt/mydatalocal/Tools/GeMoMa/GeMoMa-1.7.1.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref} a=${pathAnnotations}/${annotfile}  g=${pathGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 
fi

#########################################################################
