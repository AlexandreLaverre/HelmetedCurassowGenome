#!/bin/bash

export target=$1
export assembly=$2
export ref=$3
export source=$4
export cluster=$5
export threads=$6

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
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools
    export version=1.8
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathAllGenomes=${path}/data/genome_sequences
export pathSourceGenomes=${path}/data/genome_sequences/${source}
export pathSourceAnnotations=${path}/data/genome_annotations/${source}
export pathResults=${path}/results/genome_annotation/${target}/${assembly}/GeMoMa/${ref}
export pathScripts=${path}/scripts/gene_annotation
export pathScriptsScaffoldAssembly=${path}/scripts/scaffold_assembly

#########################################################################

if [ -e ${pathResults} ]; then
    echo "results dir already there"
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
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
    export suffix=genome_sequence
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
	    if [ -e ${pathAllGenomes}/NCBI/${target}.fa ]; then
		perl ${pathScriptsScaffoldAssembly}/cleanup.fasta.names.pl --pathInput=${pathAllGenomes}/NCBI/${target}.fa --pathOutput=${pathAllGenomes}/NCBI/${target}.clean.fa
	    else
		
		echo "looking for "${pathAllGenomes}/NCBI/${target}.fa.gz
		echo "cannot find target genome file!"
		exit
	    fi
	fi
    fi
fi

#########################################################################

export genomefile=`ls ${pathSourceGenomes} | grep ${ref}'\.' | grep fa`

#########################################################################

echo "genome file "${genomefile}

#########################################################################

if [ -e ${pathResults}/final_annotation.gff ]; then
    echo "full annotation already done for this species"
    exit
fi

#########################################################################

for annotfile in `ls ${pathSourceAnnotations}/parts | grep ${ref}'\.'`
do
    export part=`echo ${annotfile} | cut -f 2 -d '.' `

    if [ -e ${pathResults}/${part} ]; then
	echo "dir output already there"
    else
	mkdir -p ${pathResults}/${part} 
    fi

    if [ -e ${pathResults}/${part}/final_annotation.gff ]; then
	echo "already done"
    else
	echo "#!/bin/bash" > ${pathScripts}/bsub_script_gemoma

	#############################################
	
	if [ ${cluster} = "pbil" ]; then
	    echo "#SBATCH --job-name=gemoma_${ref}" >>  ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${ref}_${part}.txt" >>  ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${ref}_${part}.txt" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_gemoma
	    
	    echo "singularity exec -B ${path} -B ${pathTools} ${pathTools}/basic_ubuntu.simg java -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults}/${part} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref}_${part} a=${pathSourceAnnotations}/parts/${annotfile}  g=${pathSourceGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 m=${pathTools}/mmseqs/bin/ " >> ${pathScripts}/bsub_script_gemoma
	    
	    sbatch ${pathScripts}/bsub_script_gemoma
	fi


	#############################################
	
	if [ ${cluster} = "in2p3" ]; then
	    echo "#SBATCH --job-name=gemoma_${ref}" >>  ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${ref}_${part}.txt" >>  ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${ref}_${part}.txt" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_gemoma
	    echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_gemoma
	    
	    echo "java -Xms2G -Xmx64G -Xss1G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults}/${part} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref}_${part} a=${pathSourceAnnotations}/parts/${annotfile}  g=${pathSourceGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 m=${pathTools}/mmseqs/bin/ " >> ${pathScripts}/bsub_script_gemoma
	    
	    sbatch ${pathScripts}/bsub_script_gemoma
	fi
	
	#############################################
	
	if [ ${cluster} = "cloud" ]; then
	    ## mmseqs available in PATH
	    echo "java  -Xms2G -Xmx64G -Xss1G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults}/${part} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref}_${part} a=${pathSourceAnnotations}/parts/${annotfile}  g=${pathSourceGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 " >> ${pathScripts}/bsub_script_gemoma
	    
	    chmod a+x ${pathScripts}/bsub_script_gemoma
	    ${pathScripts}/bsub_script_gemoma
	fi
	
    fi

done

#########################################################################
