#!/bin/bash

########################################################################

export target=$1
export cluster=$2
export threads=$3
export hours=$4

#########################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbildeb" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

export assembly="Trinity"

export pathResults=${path}/results/transcriptome_assembly/${target}
export pathScripts=${path}/scripts/transcriptome_annotation

#########################################################################

export prefix=CombinedORFs_Proteins

if [ -e ${pathResults}/${prefix}.dmnd ]; then
	echo "diamond database exists"
else
    diamond makedb --ignore-warnings --in ${pathResults}/${prefix}.fa -d ${pathResults}/${prefix}
fi

#########################################################################

if [ -e ${pathResults}/${prefix}_vs_self.diamond.blastp.out ]; then
    echo "already done"
else
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_diamond

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=diamond_${target}" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --output=${pathScripts}/std_output_diamond_${target}.txt" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --error=${pathScripts}/std_error_diamond_${target}.txt" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_diamond
    fi

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=diamond_${target}" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --output=${pathScripts}/std_output_diamond_${target}.txt" >>  ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --error=${pathScripts}/std_error_diamond_${target}.txt" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_diamond
	echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_diamond
    fi

    echo "diamond blastp --threads ${threads} --db ${pathResults}/${prefix} --query ${pathResults}/${prefix}.fa --out ${pathResults}/${prefix}_vs_self.diamond.blastp.out --evalue 0.001 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps " >> ${pathScripts}/bsub_script_diamond

    if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
	sbatch ${pathScripts}/bsub_script_diamond
    fi

    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_diamond
	${pathScripts}/bsub_script_diamond
    fi
fi

#########################################################################
