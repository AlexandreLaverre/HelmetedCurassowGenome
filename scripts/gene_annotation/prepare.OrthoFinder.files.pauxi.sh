#!/bin/bash

export sp="Pauxi_pauxi"
export assembly="MEGAHIT_RAGOUT"
export genome="genome_sequence_renamed"

export cluster=$2
export annot=$3
export type=$4

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools/OrthoFinder/tools
fi

export pathProteins=${path}/data/protein_sequences
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}/${genome}.fa
export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/${annot}

mkdir ${pathResults}/OrthoFinder_${type}

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
    
    ln -s ${pathProteins}/Ensembl103/primary_transcripts/${file} ${pathResults}/OrthoFinder_${type}/${sp}.fa
done 

##########################################################################

## add Penelope pileata

ln -s ${pathProteins}/NCBI/GCA_013396635.1_ASM1339663v1_protein.faa ${pathResults}/OrthoFinder_${type}/Penelope_pileata.fa

ln -s ${pathProteins}/NCBI/GCA_013399715.1_ASM1339971v1_protein.faa ${pathResults}/OrthoFinder_${type}/Alectura_lathami.fa

ln -s ${pathProteins}/NCBI/GCA_013396415.1_ASM1339641v1_protein.faa ${pathResults}/OrthoFinder_${type}/Casuarius_casuarius.fa

##########################################################################

## add hocco

if [ ${annot} = "BRAKER_Ensembl103" ]||[ ${annot} = "BRAKER_Ensembl103_bird_species" ]; then 
    gffread -S -y ${pathResults}/braker.faa -g ${pathGenomeAssembly} ${pathResults}/braker.gtf
    
    ln -s ${pathResults}/braker.faa ${pathResults}/OrthoFinder_${type}/Pauxi_pauxi.fa
fi

##########################################################################

if [ ${annot} = "GeMoMa/combined" ]; then
    if [ ${type} = "filtered" ]; then
	if [ -e ${pathResults}/primary_transcripts/filtered_predictions_formatted.faa ]; then
	    echo "primary transcripts already done"
	else
	    python ${pathTools}/primary_transcript.py ${pathResults}/filtered_predictions_formatted.faa
	fi
	
	ln -s ${pathResults}/primary_transcripts/filtered_predictions_formatted.faa ${pathResults}/OrthoFinder_${type}/Pauxi_pauxi.fa
    fi

    if [ ${type} = "final" ]; then
	if [ -e ${pathResults}/primary_transcripts/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa ]; then
	    echo "primary transcripts already done"
	else
	    python ${pathTools}/primary_transcript.py ${pathResults}/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa
	fi
	
	ln -s ${pathResults}/primary_transcripts/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa ${pathResults}/OrthoFinder_${type}/Pauxi_pauxi.fa
    fi
fi

##########################################################################
