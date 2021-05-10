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
export pathResults=${path}/results/genome_annotation/${method}/${brakerset}

mkdir ${pathResults}/diamond_results

##########################################################################

for file in `ls ${pathProteins}/Ensembl103/primary_transcripts/ | grep all.fa`
do
    export sp=`echo ${file} | cut -f 1 -d '.'`
    
    diamond blastp --threads 24 --evalue 0.001 --max-target-seqs 20 --in ${pathResults}/braker.faa --db ${pathProteins}/Ensembl103/primary_transcripts/${sp} --out ${pathResults}/diamond_results/${sp}.diamond.blastp.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps
done 

##########################################################################

diamond blastp --threads 24 --evalue 0.001  --max-target-seqs 20  --in ${pathResults}/braker.faa --db ${pathProteins}/B10K_NCBI/Penelope_pileata --out ${pathResults}/diamond_results/Penelope_pileata.diamond.blastp.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps

##########################################################################
