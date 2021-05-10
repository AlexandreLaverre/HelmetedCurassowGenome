#!/bin/bash

export cluster=$1
export minid=$2
export minquerycov=$3
export minsubjectcov=$4

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

    if [ -e ${pathResults}/diamond_results/${sp}_minid${minid}_minquerycov${minquerycov}_minsubjectcov${minsubjectcov}.diamond.blastp.out ]; then
	echo ${sp} "already done"
    else
    	diamond blastp --threads 24 --evalue 0.001 --max-target-seqs 20 --id ${minid} --query-cover ${minquerycov} --subject-cover ${minsubjectcov} --query ${pathResults}/braker.faa --db ${pathProteins}/Ensembl103/primary_transcripts/${sp} --out ${pathResults}/diamond_results/${sp}_minid${minid}_minquerycov${minquerycov}_minsubjectcov${minsubjectcov}.diamond.blastp.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps
    fi
done 

##########################################################################

if [ -e ${pathResults}/diamond_results/Penelope_pileata_minid${minid}_minquerycov${minquerycov}_minsubjectcov${minsubjectcov}.diamond.blastp.out ]; then
    echo "Penelope pileata already done"
else
    diamond blastp --threads 24 --evalue 0.001  --max-target-seqs 20  --id ${minid} --query-cover ${minquerycov} --subject-cover ${minsubjectcov}  --query ${pathResults}/braker.faa --db ${pathProteins}/B10K_NCBI/Penelope_pileata --out ${pathResults}/diamond_results/Penelope_pileata_minid${minid}_minquerycov${minquerycov}_minsubjectcov${minsubjectcov}.diamond.blastp.out --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps
fi

##########################################################################
