#!/bin/bash

export cluster=$1

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathCDS=${path}/data/coding_sequences
export pathAnnot=${path}/results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined
export pathResults=${path}/results/coding_gene_evolution
export pathScripts=${path}/scripts/coding_gene_evolution

##########################################################################

for file in `ls ${pathResults}/CDS | grep unaln `
do
    export prefix=`basename ${file} .unaln.fa`
    
    prank -codon -f=fasta -t=${pathResults}/CDS/${prefix}.tree -d=${pathResults}/CDS/${file} -o=${pathResults}/CDS/${prefix}.aln.fa

    prank -convert -f=phylips -d=${pathResults}/CDS/${prefix}.aln.fa -o=${pathResults}/CDS/${prefix}.aln.phy
    
done

##########################################################################
