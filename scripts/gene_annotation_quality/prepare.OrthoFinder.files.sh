#!/bin/bash

export cluster=$1
export annot=$2

export assembly="MEGAHIT_RAGOUT"
export genome="genome_sequence_renamed"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathProteins=${path}/data/protein_sequences
export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}/${genome}.fa
export pathResults=${path}/results/genome_annotation/${assembly}/${annot}

mkdir ${pathResults}/OrthoFinder

##########################################################################

for file in `ls ${pathProteins}/Ensembl103 | grep all.fa`
do
    if [ -e ${pathProteins}/Ensembl103/primary_transcripts/${file} ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathProteins}/Ensembl103/${file}
    fi
done

##########################################################################

for file in `ls ${pathProteins}/Ensembl103/primary_transcripts/ | grep all.fa`
do
    export sp=`echo ${file} | cut -f 1 -d '.'`
    
    ln -s ${pathProteins}/Ensembl103/primary_transcripts/${file} ${pathResults}/OrthoFinder/${sp}.fa
done 

##########################################################################

## add Penelope pileata

ln -s ${pathProteins}/B10K_NCBI/GCA_013396635.1_ASM1339663v1_protein.faa ${pathResults}/OrthoFinder/Penelope_pileata.fa

##########################################################################

## add hocco

if [ ${annot} = "BRAKER_Ensembl103_multithread" ]||[ ${annot} = "BRAKER_Ensembl103_singlethread" ]||[ ${annot} = "BRAKER_Ensembl103_bird_species" ]; then 
    gffread -S -y ${pathResults}/braker.faa -g ${pathGenomeAssembly} ${pathResults}/braker.gtf
    
    ln -s ${pathResults}/braker.faa ${pathResults}/OrthoFinder/Pauxi_pauxi.fa
fi

##########################################################################

if [ ${annot} = "GeMoMa/combined" ]; then 
    if [ -e ${pathResults}/primary_transcripts/filtered_predictions_formatted.faa ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathResults}/filtered_predictions_formatted.faa
    fi
    
    ln -s ${pathResults}/primary_transcripts/filtered_predictions_formatted.faa ${pathResults}/OrthoFinder/Pauxi_pauxi.fa

fi

##########################################################################
