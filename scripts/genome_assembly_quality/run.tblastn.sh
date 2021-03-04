#!/bin/bash

########################################################################

export sp="Chicken"

export method=$1	# ie : SOAPdenovo ABYSS MEGAHIT IDBA Discovar
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

export pathProteinSequences=${path}/data/protein_sequences/${sp}
export pathResults=${path}/results/genome_assembly_tests/${method}/${genomeprefix} 
export pathScripts=${path}/scripts/genome_assembly

export ensrelease=98

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

## first make blastdb

if [ -e ${pathResults}/${suffix}.nhr ]; then
    echo "blast database already done"
else
    makeblastdb -dbtype nucl -in ${pathAssembly} -out ${pathResults}/${suffix}
fi

#########################################################################

if [ -e ${pathResults}/AllPeptides${ensrelease}_vs_${suffix}.tblastn.out ]; then
    echo "already done"
else
    tblastn -num_threads ${threads} -query ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa -db ${pathResults}/${suffix} -out ${pathResults}/AllPeptides${ensrelease}_vs_${suffix}.tblastn.out -evalue 0.001 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps"
fi

#########################################################################
