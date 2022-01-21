#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export refsp=$3
export cluster=$4
export nthreads=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathProteinSequences=${path}/data/protein_sequences/${refsp}
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathResults=${path}/results/genome_assembly_quality/${sp}/${assembly}
export pathScripts=${path}/scripts/genome_assembly

export ensrelease=103

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence.fa
    export suffix=genome_sequence
fi

#########################################################################

if [ ${nthreads} = "parts" ]; then
    for i in {0..100};
    do
	if [ -e ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa ]; then
	    export last=`tail -n 1 ${pathResults}/tblastn_parts/${refsp}_AllPeptides${ensrelease}_vs_${suffix}_part${i}.tblastn.out | cut -f 1`
	    export tot=`grep -c ">" ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa `
	    export index=`grep ">" ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa | grep -n ${last} | cut -f 1 -d ':'`
	    
	    export ratio=$(($index * 100 / $tot))
    
	    echo "part "${i}" index "${index}" out of "${tot}" "${ratio}"% done";
	fi

    done
else
    export last=`tail -n 1 ${pathResults}/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out | cut -f 1`
    export tot=`grep -c ">" ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa `
    export index=`grep ">" ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa | grep -n ${last} | cut -f 1 -d ':'`
    
    export ratio=$(($index * 100 / $tot))
    
    echo "index "${index}" out of "${tot}" "${ratio}"% done";
fi
#########################################################################
