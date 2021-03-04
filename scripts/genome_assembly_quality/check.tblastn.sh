#!/bin/bash

########################################################################

export sp="Chicken"

export method=$1        # ie : SOAPdenovo ABYSS MEGAHIT IDBA Discovar
export kmer=$2          # ie : 50, 79...
export genomeprefix=$3  # ie : HelmetedCurassow ...
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

export pathProteinSequences=${path}/data/protein_sequences/${sp}
export pathResults=${path}/results/genome_assembly_tests/${method}/${genomeprefix}
export pathScripts=${path}/scripts/genome_assembly

export ensrelease=98

#########################################################################

if [ ${method} = "SOAPdenovo" ]; then
    export suffix=kmer${kmer}.scafSeq
fi

if [ ${method} = "ABYSS" ]; then
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

export last=`tail -n 1 ${pathResults}/AllPeptides${ensrelease}_vs_${suffix}.tblastn.out | cut -f 1`
export tot=`grep -c ">" ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa `
export index=`grep ">" ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa | grep -n ${last} | cut -f 1 -d ':'`

export ratio=$(($index * 100 / $tot))

echo "index "${index}" out of "${tot}" "${ratio}"% done";

#########################################################################
