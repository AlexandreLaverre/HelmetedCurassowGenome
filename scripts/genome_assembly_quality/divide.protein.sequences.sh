#!/bin/bash

export sp=$1
export nbparts=$2
export cluster=$3

#####################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
else
    if [ ${cluster} = "cloud" ]; then
	export path=/mnt/mydatalocal/HelmetedCurassowGenome/
    else
	echo "unknown cluster"
	exit
    fi
fi

export pathProteinSequences=${path}/data/protein_sequences/${sp}
export pathScripts=${path}/scripts/genome_assembly_quality

export release=103

###################################################################################

if [ -e ${pathProteinSequences}/fasta_parts ]; then
    echo "path output already there"
else
   mkdir ${pathProteinSequences}/fasta_parts
fi

###################################################################################

perl ${pathScripts}/divide.sequences.pl --pathFastaInput=${pathProteinSequences}/AllPeptides_Ensembl${release}.fa --nbParts=${nbparts} --prefixOutput=${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${release}

###################################################################################

