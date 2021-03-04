#!/bin/bash

########################################################################

export sp="Chicken"

export method=$1	# ie : SOAPdenovo ABYSS MEGAHIT IDBA Discovar
export kmer=$2		# ie : 50, 79...
export genomeprefix=$3  # ie : HelmetedCurassow ...
export user=$4		# ie : necsulea or alaverre
export cluster=$5	# ie : pbil or cloud
export threads=$6	# ie : number of threads

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

if [ ${method} = "Discovar" ]; then
    export pathAssembly=${pathResults}/a.final/a.fasta
    export suffix=final.assembly
fi

if [ ${method} = "IDBA" ]; then
    export pathAssembly=${pathResults}/scaffold.fa
    export suffix=scaffolds
fi

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathResults}/final.contigs.fa
    export suffix=final.contigs
fi

if [ ${method} = "Velvet" ]; then
    export pathAssembly=${pathResults}/_${kmer}/contigs.fa
    export suffix=kmer${kmer}.contigs
fi

#########################################################################

pblat -threads=${threads} -t=dna -q=dna -out=blast8 ${pathGenome}/${genomeprefix}.fa ${pathAssembly} ${pathResults}/${suffix}.blat_vs_${genomeprefix}.out

#########################################################################
