#!/bin/bash

export target=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
    export pathTools=/beegfs/home/${USER}/Tools
    export version=1.8
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
    export pathTools=/sps/biometr/necsulea/Tools
    export version=1.8
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools
    export version=1.8
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${target}/${assembly}
export pathAllGenomes=${path}/data/genome_sequences
export pathRNASeq=${path}/results/RNASeq_alignments/${target}/${assembly}
export pathResults=${path}/results/genome_annotation/${target}/${assembly}/GeMoMa
export pathScripts=${path}/scripts/gene_annotation
export pathScriptsScaffoldAssembly=${path}/scripts/scaffold_assembly

## using mmseqs2, installed using apt on Ubuntu 20.04
## version 9-d36de+ds-4 amd64

#########################################################################

if [ -e ${pathResults} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    if [ -e ${pathGenomeAssembly}/genome_sequence_renamed_sm.fa ]; then
	export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
    else
	if [ -e ${pathGenomeAssembly}/genome_sequence_renamed.fa ]; then
	    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
	else
	    echo "cannot find genome sequence"
	    exit
	fi
    fi
fi

#########################################################################

if [ ${assembly} = "NCBI" ]; then
    export pathAssembly=${pathAllGenomes}/NCBI/${target}.clean.fa

    if [ -e ${pathAssembly} ]; then
	echo "genome file already there"
    else
	if [ -e ${pathAllGenomes}/NCBI/${target}.fa.gz ]; then
	    perl ${pathScriptsScaffoldAssembly}/cleanup.fasta.names.pl --pathInput=${pathAllGenomes}/NCBI/${target}.fa.gz --pathOutput=${pathAllGenomes}/NCBI/${target}.clean.fa
	else
	    echo "looking for "${pathAllGenomes}/NCBI/${target}.fa.gz
	    echo "cannot find target genome file!"
	    exit
	fi
    fi
fi

#########################################################################

if [ -e ${pathResults}/final_annotation.gff ]; then
    echo "already done"
else
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_gemoma

    #############################################

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=gemoma_rna_seq_${target}" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${target}.txt" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${target}.txt" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_gemoma

	echo "singularity exec -B ${path} -B ${pathTools} ${pathTools}/basic_ubuntu.simg java -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI ERE s=FR_FIRST_STRAND m=${pathRNASeq}/accepted_hits_all_samples.bam outdir=${pathResults}" >> ${pathScripts}/bsub_script_gemoma

	sbatch ${pathScripts}/bsub_script_gemoma
    fi


    #############################################

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=gemoma_${target}" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${target}.txt" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${target}.txt" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_gemoma

	echo "java -Xms2G -Xmx64G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI ERE s=FR_FIRST_STRAND m=${pathRNASeq}/accepted_hits_all_samples.bam outdir=${pathResults}" >> ${pathScripts}/bsub_script_gemoma

	sbatch ${pathScripts}/bsub_script_gemoma
    fi

    #############################################

    if [ ${cluster} = "cloud" ]; then
	## mmseqs available in PATH
	echo "java  -Xms2G -Xmx64G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI ERE s=FR_FIRST_STRAND m=${pathRNASeq}/accepted_hits_all_samples.bam outdir=${pathResults}" >> ${pathScripts}/bsub_script_gemoma

	chmod a+x ${pathScripts}/bsub_script_gemoma
	${pathScripts}/bsub_script_gemoma
    fi

fi

#########################################################################
