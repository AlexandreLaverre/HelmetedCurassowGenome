#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export refsp=$3
export source=$4
export cluster=$5
export threads=$6
export hours=$7

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
export pathResults=${path}/results/genome_annotation/${target}/${assembly}/GeMoMa/combined

## diamond 2.0.15

#########################################################################

export protfile=`ls ${pathProteinSequences} | grep "${refsp}\." | grep fa`
echo ${protfile}

if [ -e ${pathProteinSequences}/${protfile} ]; then
    echo "using as an input "${pathProteinSequences}/${protfile}

    export prefix=`basename ${protfile} .fa`

    if [ -e ${pathProteinSequences}/${prefix}.dmnd ]; then
	echo "diamond database exists"
    else
	diamond makedb --in ${pathProteinSequences}/${protfile} -d ${pathProteinSequences}/${prefix}
    fi
else
    echo "cannot find protein sequence file for "${refsp}
    exit
fi

#########################################################################

if [ -e ${pathResults}/diamond_results ]; then
    echo "output dir already there"
else
    mkdir -p ${pathResults}/diamond_results
fi

#########################################################################

if [ -e ${pathResults}/GeMoMa_vs_${refsp}.diamond.blastp.out ]; then
    echo "already done"
else
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_diamond

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=diamond_${refsp}" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --output=${pathScripts}/std_output_diamond_${refsp}.txt" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --error=${pathScripts}/std_error_diamond_${refsp}.txt" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_diamond
    fi

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=diamond_${refsp}" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --output=${pathScripts}/std_output_diamond_${refsp}.txt" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --error=${pathScripts}/std_error_diamond_${refsp}.txt" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_diamond
    fi

    echo "diamond blastp --db ${pathProteinSequences}/${prefix} --query ${pathResults}/filtered_predictions_formatted.faa --out ${pathResults}/diamond_results/GeMoMa_vs_${refsp}.diamond.blastp.out --evalue 0.001 --outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps\"" >> ${pathScripts}/bsub_script_diamond

    if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
	sbatch ${pathScripts}/bsub_script_diamond
    fi

    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_diamond
	${pathScripts}/bsub_script_diamond
    fi
fi

#########################################################################
