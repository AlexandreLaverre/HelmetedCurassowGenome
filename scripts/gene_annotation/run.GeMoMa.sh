#!/bin/bash

export assembly=$1
export ref=$2
export source=$3
export cluster=$4
export threads=$5
export constraint=$6

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
export pathGenomes=${path}/data/genome_sequences/${source}
export pathAnnotations=${path}/data/genome_annotations/${source}
export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/${ref}
export pathScripts=${path}/scripts/gene_annotation

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

if [ ${source} = "Ensembl" ]; then
    export genomefile=`ls ${pathGenomes} | grep ${ref}'\.' | grep fa`
    export annotfile=`ls ${pathAnnotations} | grep ${ref}'\.' | grep gff`
fi

if [ ${source} = "NCBI" ]; then
    export genomefile=${ref}.fa
    export annotfile=${ref}_withstopcodons.gff
fi

#########################################################################

echo "genome file "${genomefile}
echo "annot file "${annotfile}

#########################################################################

if [ -e ${pathResults}/final_annotation.gff ]; then
    echo "already done"
else
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_gemoma

    
    #############################################

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=gemoma_${ref}" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${ref}.txt" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${ref}.txt" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --constraint=${constraint}">> ${pathScripts}/bsub_script_gemoma

	echo "singularity exec -B ${path} -B ${pathTools} ${pathTools}/basic_ubuntu.simg java -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref} a=${pathAnnotations}/${annotfile}  g=${pathGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 m=${pathTools}/mmseqs/bin/ " >> ${pathScripts}/bsub_script_gemoma

	sbatch ${pathScripts}/bsub_script_gemoma
    fi


    #############################################

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=gemoma_${ref}" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${ref}.txt" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${ref}.txt" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_gemoma
	
	echo "java -Xms2G -Xmx64G  -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref} a=${pathAnnotations}/${annotfile}  g=${pathGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 m=${pathTools}/mmseqs/bin/ " >> ${pathScripts}/bsub_script_gemoma

	sbatch ${pathScripts}/bsub_script_gemoma
    fi

    
    #############################################
    
    if [ ${cluster} = "cloud" ]; then
	## mmseqs available in PATH
	echo "java  -Xms2G -Xmx64G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathAssembly} i=${ref} a=${pathAnnotations}/${annotfile}  g=${pathGenomes}/${genomefile} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 " >> ${pathScripts}/bsub_script_gemoma
	
	chmod a+x ${pathScripts}/bsub_script_gemoma
	${pathScripts}/bsub_script_gemoma
    fi
    
fi

#########################################################################
