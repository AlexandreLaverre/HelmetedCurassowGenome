#!/bin/bash

########################################################################

export method=$1
export refsp=$2
export cluster=$3 

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathProteinSequences=${path}/data/protein_sequences/${refsp}
export pathGenomeAssembly=${path}/results/genome_assembly/${method}
export pathResults=${path}/results/genome_assembly_quality/${method}
export pathScripts=${path}/scripts/genome_assembly

export ensrelease=103

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

perl ${pathScripts}/analyze.tblastn.statistics.pl --pathProteins=${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa pathTBlastNResults=${pathResults}/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out --minPCIdentity=${minPCIdentity} --maxEValue=${maxEValue} --pathOutput=${pathResults}/AlignmentStatistics_${refsp}_AllPeptides${ensrelease}_vs_${suffix}.txt

#########################################################################
