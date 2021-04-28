#!/bin/bash

########################################################################

export method=$1
export refsp=$2
export cluster=$3 
export threads=$4

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

## first make blastdb

if [ -e ${pathResults}/${suffix}.nhr ]; then
    echo "blast database already done"
else
    makeblastdb -dbtype nucl -in ${pathAssembly} -out ${pathResults}/${suffix}
fi

#########################################################################

if [ ${threads} = "parts" ]; then

    if [ -e ${pathResults}/tblastn_parts ]; then
	echo "output parts already there"
    else
	mkdir ${pathResults}/tblastn_parts
    fi
    
    for i in {0..100}
    do
	if [ -e ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa ]; then
	    tblastn -num_threads 1 -query ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa -db ${pathResults}/${suffix} -out ${pathResults}/tblastn_parts/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out -evalue 0.001 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps" & 
	fi
    done
else
    
    if [ -e ${pathResults}/AllPeptides${ensrelease}_vs_${suffix}.tblastn.out ]; then
	echo "already done"
    else
	tblastn -num_threads ${threads} -query ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa -db ${pathResults}/${suffix} -out ${pathResults}/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out -evalue 0.001 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps"
    fi
fi

#########################################################################
