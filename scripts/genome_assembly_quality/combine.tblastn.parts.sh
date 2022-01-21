#!/bin/bash

########################################################################

export sp=$1
export method=$2
export refsp=$3
export cluster=$4

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathProteinSequences=${path}/data/protein_sequences/${refsp}
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${method}
export pathResults=${path}/results/genome_assembly_quality/${sp}/${method}
export pathScripts=${path}/scripts/genome_assembly

export ensrelease=103

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence.fa
    export suffix=genome_sequence
fi

#########################################################################

for i in {0..100}
do
    if [ -e ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa ]; then
	
	cat ${pathResults}/tblastn_parts/${refsp}_AllPeptides${ensrelease}_vs_${suffix}_part${i}.tblastn.out >> ${pathResults}/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out 
	
    fi
    
done

#########################################################################
