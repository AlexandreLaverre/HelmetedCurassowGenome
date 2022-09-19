#!/bin/bash

export cluster=$1

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
    export pathTools=/beegfs/home/necsulea/OrthoFinder/tools
fi

export pathEnsemblProteins=${path}/data/protein_sequences/Ensembl103
export pathGeMoMa=${path}/results/genome_annotation
export pathResults=${path}/results/gene_families
export pathScripts=${path}/scripts/gene_families

##########################################################################

if [ -e ${pathResults} ]; then
    echo "output path already there"
else
    mkdir -p ${pathResults}
fi

##########################################################################

## Ensembl species

for sp in `grep Ensembl ${pathScripts}/species_list.txt | cut -f 1`
do
    export file=`ls ${pathEnsemblProteins}/| grep all.fa | grep ${sp}'\.'`

    if [ -e ${pathEnsemblProteins}/primary_transcripts/${file} ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathEnsemblProteins}/${file}
    fi

    ln -s ${pathEnsemblProteins}/primary_transcripts/${file} ${pathResults}/${sp}.fa
done

##########################################################################

## NCBI species

for sp in `grep NCBI ${pathScripts}/species_list.txt | cut -f 1`
do
    if [ -e ${pathGeMoMa}/${sp}/NCBI/GeMoMa/combined/primary_transcripts/combined_annotations_NCBI_GeMoMa.faa ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathGeMoMa}/${sp}/NCBI/GeMoMa/combined/combined_annotations_NCBI_GeMoMa.faa
    fi
    
    ln -s ${pathGeMoMa}/${sp}/NCBI/GeMoMa/combined/primary_transcripts/combined_annotations_NCBI_GeMoMa.faa ${pathResults}/${sp}.fa
done

##########################################################################

## Basiliscus and Pauxi

for sp in Pauxi_pauxi Basiliscus_vittatus

if [ -e ${pathGeMoMa}/${sp}/MEGAHIT_RAGOUT/GeMoMa/combined/primary_transcripts/filtered_GeMoMa_annotations.faa ]; then
    echo "primary transcripts already done"
else
    python ${pathTools}/primary_transcript.py ${pathGeMoMa}/${sp}/MEGAHIT_RAGOUT/GeMoMa/combined/filtered_GeMoMa_annotations.faa
fi

ln -s ${pathGeMoMa}/${sp}/MEGAHIT_RAGOUT/GeMoMa/combined/primary_transcripts/filtered_GeMoMa_annotations.faa ${pathResults}/${sp}.fa

##########################################################################
