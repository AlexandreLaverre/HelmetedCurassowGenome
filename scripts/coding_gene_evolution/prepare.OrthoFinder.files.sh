#!/bin/bash

export cluster=$1

export assembly="MEGAHIT_RAGOUT"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathProteins=${path}/data/protein_sequences
export pathResults=${path}/results/coding_sequence_evolution
export pathAnnot=${path}/results/genome_annotation/GeMoMa/combined
export pathScripts=${path}/scripts/coding_sequence_evolution

##########################################################################

for sp in `grep Ensembl ${pathScripts}/species_list | cut -f 1`
do
    export file=`ls ${pathProteins}/Ensembl103 | grep all.fa | grep ${sp}'\.'`

    if [ -e ${pathProteins}/Ensembl103/primary_transcripts/${file} ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathProteins}/Ensembl103/${file}
    fi

    ln -s ${pathProteins}/Ensembl103/primary_transcripts/${file} ${pathResults}/${sp}.fa
done

##########################################################################

## add NCBI species

ln -s ${pathProteins}/NCBI/GCA_013396635.1_ASM1339663v1_protein.faa ${pathResults}/Penelope_pileata.fa

ln -s ${pathProteins}/NCBI/GCA_013399715.1_ASM1339971v1_protein.faa ${pathResults}/Alectura_lathami.fa

ln -s ${pathProteins}/NCBI/GCA_013396415.1_ASM1339641v1_protein.faa ${pathResults}/Casuarius_casuarius.fa

ln -s ${pathProteins}/NCBI/GCA_013399115.1_ASM1339911v1_protein.faa ${pathResults}/Anseranas_semipalmata.fa

ln -s ${pathProteins}/NCBI/GCF_000709895.1_ASM70989v1_protein.faa ${pathResults}/Balearica_regulorum.fa

ln -s ${pathProteins}/NCBI/GCA_013390085.1_ASM1339008v1_protein.faa ${pathResults}/Grus_americana.fa

ln -s ${pathProteins}/NCBI/GCA_013398885.1_ASM1339888v1_protein.faa ${pathResults}/Bucorvus_abyssinicus.fa

ln -s ${pathProteins}/NCBI/GCA_000710305.1_ASM71030v1_protein.faa ${pathResults}/Buceros_rhinoceros.fa

ln -s ${pathProteins}/NCBI/GCA_013397515.1_ASM1339751v1_protein.faa ${pathResults}/Upupa_epops.fa

ln -s ${pathProteins}/NCBI/GCA_013400115.1_ASM1340011v1_protein.faa ${pathResults}/Rhinopomastus_cyanomelas.fa

##########################################################################

## add hocco

if [ -e ${pathAnnot}/primary_transcripts/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa ]; then
    echo "primary transcripts already done"
else
    python ${pathTools}/primary_transcript.py ${pathAnnot}/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa
fi

ln -s ${pathAnnot}/primary_transcripts/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa ${pathResults}/Pauxi_pauxi.fa

##########################################################################
