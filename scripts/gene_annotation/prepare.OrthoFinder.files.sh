#!/bin/bash

export cluster=$1
export method="MEGAHIT_RAGOUT"
export brakerset="BRAKER_Ensembl103_multithread"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathProteins=${path}/data/protein_sequences
export pathResults=${path}/results/genome_annotation/${method}/${brakerset}

mkdir ${pathResults}/OrthoFinder

##########################################################################

# for file in `ls ${pathProteins}/Ensembl103 | grep all.fa`
# do
#     python ${pathTools}/primary_transcript.py ${pathProteins}/Ensembl103/${file}
# done

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
