#!/bin/bash

########################################################################

export sp="Chicken"

export method=$1
export kmer=$2
export genomeprefix=$3
export user=$4          # ie : necsulea or alaverre
export cluster=$5       # ie : pbil or cloud

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${user}/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    if [ ${user} = "necsulea" ]; then
	export path=/mnt/IPLOSS; else
	export path=/mnt/
    fi
fi

export pathGenome=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/genome_assembly_tests/${method}/${genomeprefix}
export pathScripts=${path}/scripts/genome_assembly

#########################################################################

if [ ${method} = "SOAPdenovo" ]; then
    export pathAssembly=${pathResults}/kmer${kmer}.scafSeq
    export suffix=kmer${kmer}.scafSeq
fi


if [ ${method} = "ABYSS" ]; then
    export pathAssembly=${pathResults}/kmer${kmer}-scaffolds.fa
    export suffix=kmer${kmer}-scaffolds
fi


#########################################################################

blat -t=dna -q=dna -out=blast8 ${pathGenome}/${genomeprefix}.fa ${pathAssembly} ${pathResults}/${suffix}.blat_vs_${genomeprefix}.out

#########################################################################
