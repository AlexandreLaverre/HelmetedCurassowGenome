#!/bin/bash

export cluster=$1
export method="MEGAHIT_RAGOUT"
export brakerset="BRAKER_Ensembl103_multithread"
export genome="genome_sequence_renamed"

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi

export pathProteins=${path}/data/protein_sequences

##########################################################################

for file in `ls ${pathProteins}/Ensembl103/primary_transcripts/ | grep all.fa`
do
    export sp=`echo ${file} | cut -f 1 -d '.'`
    
    diamond makedb --in ${pathProteins}/Ensembl103/primary_transcripts/${file} --db ${pathProteins}/Ensembl103/primary_transcripts/${sp}
done 

##########################################################################

diamond makedb --in ${pathProteins}/B10K_NCBI/GCA_013396635.1_ASM1339663v1_protein.faa ${pathProteins}/B10K_NCBI/Penelope_pileata 

##########################################################################
