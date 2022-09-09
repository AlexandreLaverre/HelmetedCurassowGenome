#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export refsp=$3
export source=$4
export cluster=$5
export threads=$6

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

export pathProteinSequences=${path}/data/protein_sequences/${source}/primary_transcripts
export pathTranscriptomeAssembly=${path}/results/transcriptome_assembly/${target}
export pathResults=${path}/results/transcriptome_assembly/${target}/tblastn_results
export pathScripts=${path}/scripts/transcriptome_annotation

## ncbi-blast-2.8.1+ on pbil
## ncbi-blast 2.12.0+ on cloud
## ncbi-blast 2.13.0+ on in2p3

#########################################################################

if [ -e ${pathResults} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}
fi

#########################################################################

if [ ${assembly} = "Trinity" ]; then
    export pathAssembly=${pathTranscriptomeAssembly}/Trinity.fasta
    export suffix=Trinity
fi

#########################################################################

export protfile=`ls ${pathProteinSequences} | grep ${refsp} | grep fa`

if [ -e ${pathProteinSequences}/${protfile} ]; then
    echo "using as an input "${pathProteinSequences}/${protfile}
else
    echo "cannot find protein sequence file for "${refsp}
    exit
fi

#########################################################################

## first make blastdb if not already there

if [ -e ${pathTranscriptomeAssembly}/${suffix}.nhr ]; then
    echo "blast database already done"
else
    echo "cannot find path" ${pathTranscriptomeAssembly}/${suffix}.nhr
    makeblastdb -dbtype nucl -in ${pathAssembly} -out ${pathTranscriptomeAssembly}/${suffix}
fi

#########################################################################

if [ -e ${pathResults}/${refsp}_vs_${suffix}.tblastn.out ]; then
    echo "already done"
else
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_tblastn

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=tblastn_${refsp}" >>  ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --output=${pathScripts}/std_output_tblastn_${refsp}.txt" >>  ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --error=${pathScripts}/std_error_tblastn_${refsp}.txt" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_tblastn
    fi

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=tblastn_${refsp}" >>  ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --output=${pathScripts}/std_output_tblastn_${refsp}.txt" >>  ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --error=${pathScripts}/std_error_tblastn_${refsp}.txt" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_tblastn
	echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_tblastn
    fi

    echo "tblastn -num_threads ${threads} -query ${pathProteinSequences}/${protfile} -db ${pathTranscriptomeAssembly}/${suffix} -out ${pathResults}/${refsp}_vs_${suffix}.tblastn.out -evalue 0.001 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps\"">> ${pathScripts}/bsub_script_tblastn

    if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
	sbatch ${pathScripts}/bsub_script_tblastn
    fi

    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_tblastn
	${pathScripts}/bsub_script_tblastn
    fi
fi

#########################################################################
