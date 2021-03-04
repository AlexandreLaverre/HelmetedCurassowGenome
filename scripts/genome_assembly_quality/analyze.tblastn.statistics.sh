#!/bin/bash

########################################################################

export sp="Chicken"

export method=$1	# ie : SOAPdenovo ABYSS MEGAHIT IDBA Discovar
export kmer=$2		# ie : 50, 79...
export genomeprefix=$3  # ie : HelmetedCurassow ...
export minPCIdentity=$4
export maxEValue=$5
export user=$6		# ie : necsulea or alaverre
export cluster=$7	# ie : pbil or cloud

export ensrelease=98

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

export pathProteinSequences=${path}/data/protein_sequences/${sp}
export pathResults=${path}/results/genome_assembly/${method}/${genomeprefix} #_tests
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

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathResults}/final.contigs.fa
    export suffix=final.contigs.fa
fi

if [ ${method} = "IDBA" ]; then
    export pathAssembly=${pathResults}/scaffold.fa
    export suffix=scaffold.fa
fi

if [ ${method} = "Discovar" ]; then
    export pathAssembly=${pathResults}/a.final/a.fasta
    export suffix=a.fasta
fi

#########################################################################

perl ${pathScripts}/analyze.tblastn.statistics.pl ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa ${pathResults}/AllPeptides${ensrelease}_vs_${suffix}.tblastn.out ${minPCIdentity} ${maxEValue} ${pathResults}/

#########################################################################
