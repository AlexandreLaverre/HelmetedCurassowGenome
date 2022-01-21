#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export refsp=$3
export cluster=$4 
export threads=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathProteinSequences=${path}/data/protein_sequences/${refsp}
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathResults=${path}/results/genome_assembly_quality/${sp}/${assembly}
export pathScripts=${path}/scripts/genome_assembly_quality

export ensrelease=103

## ncbi-blast-2.8.1+ on pbil

#########################################################################

if [ -e ${pathResults} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}
fi

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

	    echo "#!/bin/bash" > ${pathScripts}/bsub_script_tblastn
	    
	    if [ ${cluster} = "pbil" ]; then
		echo "#SBATCH --job-name=tblastn_${ref}_part${i}" >>  ${pathScripts}/bsub_script_tblastn
		echo "#SBATCH --output=${pathScripts}/std_output_tblastn_${ref}.txt" >>  ${pathScripts}/bsub_script_tblastn
		echo "#SBATCH --error=${pathScripts}/std_error_tblastn_${ref}.txt" >> ${pathScripts}/bsub_script_tblastn
		echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_tblastn
		echo "#SBATCH --mem=5G" >> ${pathScripts}/bsub_script_tblastn
		echo "#SBATCH --cpus-per-task=2" >> ${pathScripts}/bsub_script_tblastn
		echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_tblastn
	    fi
	    
	    echo "tblastn -num_threads 2 -query ${pathProteinSequences}/fasta_parts/AllPeptides_Ensembl${ensrelease}_part${i}.fa -db ${pathResults}/${suffix} -out ${pathResults}/tblastn_parts/${refsp}_AllPeptides${ensrelease}_vs_${suffix}_part${i}.tblastn.out -evalue 0.001 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps\"" >> ${pathScripts}/bsub_script_tblastn

	    if [ ${cluster} = "pbil" ]; then
		sbatch ${pathScripts}/bsub_script_tblastn
	    fi

	    if [ ${cluster} = "cloud" ]; then
		chmod a+x ${pathScripts}/bsub_script_tblastn
		${pathScripts}/bsub_script_tblastn
	    fi
	    
	fi
    done
else
    
    if [ -e ${pathResults}/AllPeptides${ensrelease}_vs_${suffix}.tblastn.out ]; then
	echo "already done"
    else
	echo "#!/bin/bash" > ${pathScripts}/bsub_script_tblastn
	
	 if [ ${cluster} = "pbil" ]; then
	     echo "#SBATCH --job-name=tblastn_${ref}" >  ${pathScripts}/bsub_script_tblastn
	     echo "#SBATCH --output=${pathScripts}/std_output_tblastn_${ref}.txt" >>  ${pathScripts}/bsub_script_tblastn
	     echo "#SBATCH --error=${pathScripts}/std_error_tblastn_${ref}.txt" >> ${pathScripts}/bsub_script_tblastn
	     echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_tblastn
	     echo "#SBATCH --mem=5G" >> ${pathScripts}/bsub_script_tblastn
	     echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_tblastn
	     echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_tblastn
	 fi
	 
	
	 echo "tblastn -num_threads ${threads} -query ${pathProteinSequences}/AllPeptides_Ensembl${ensrelease}.fa -db ${pathResults}/${suffix} -out ${pathResults}/${refsp}_AllPeptides${ensrelease}_vs_${suffix}.tblastn.out -evalue 0.001 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps\"">> ${pathScripts}/bsub_script_tblastn

	 if [ ${cluster} = "pbil" ]; then
	     sbatch ${pathScripts}/bsub_script_tblastn
	 fi
	 
	 if [ ${cluster} = "cloud" ]; then
	     chmod a+x ${pathScripts}/bsub_script_tblastn
	     ${pathScripts}/bsub_script_tblastn
	 fi
    fi
fi

#########################################################################
