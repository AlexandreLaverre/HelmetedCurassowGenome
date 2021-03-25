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

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence.fa
    export suffix=genome_sequence
fi

#########################################################################

export last=`tail -n 1 ${pathResults}/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out | cut -f 1`
export tot=`grep -c ">" ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa `
export index=`grep ">" ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa | grep -n ${last} | cut -f 1 -d ':'`

export ratio=$(($index * 100 / $tot))

echo "index "${index}" out of "${tot}" "${ratio}"% done";

#########################################################################
