#!/bin/bash

########################################################################

export sp=$1
export method=$2	# MEGAHIT or MEGAHIT_RAGOUT
export cluster=$3	# pbil or cloud
export threads=$4

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
    export samtools="singularity exec -B /beegfs/data/necsulea /beegfs/home/necsulea/Tools/samtools_1.3.1.sif samtools"
fi

export pathWGS=${path}/data/WGS/${sp}
export pathResults=${path}/results/genome_assembly/${sp}/${method}
export pathScripts=${path}/scripts/genome_assembly_quality

#########################################################################

if [ ${method} = "MEGAHIT" ]; then
    export pathAssembly=${pathResults}/final.contigs.fa
    export prefix=final.contigs
fi

if [ ${method} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathResults}/genome_sequence.fa
    export prefix=genome_sequence
fi

#########################################################################

if [ -e ${pathResults}/${prefix}.1.bt2 ]; then
    echo "bowtie index already done"
else
    bowtie2-build --threads ${threads} ${pathAssembly} ${pathResults}/${prefix}
fi

#########################################################################

mkdir -p ${pathResults}/remapped_WGS_reads

#########################################################################

if [ ${sp} = "Basiliscus_vittatus" ]; then
    for sample in Liver_Female Liver_Male
    do
	export pathR1=${pathWGS}/${sample}_R1_trimmed.fastq.gz
	export pathR2=${pathWGS}/${sample}_R2_trimmed.fastq.gz

	echo "#!/bin/bash" > ${pathScripts}/bsub_script_bowtie2
	
	 if [ ${cluster} = "pbil" ]; then
	     echo "#SBATCH --job-name=bowtie2_${sample}" >  ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --output=${pathScripts}/std_output_bowtie2_${sample}.txt" >>  ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --error=${pathScripts}/std_error_bowtie2_${sample}.txt" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --mem=5G" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_bowtie2
	 fi

	 echo "bowtie2 -x ${pathResults}/${prefix} -1 ${pathR1} -2 ${pathR2} -S ${pathResults}/remapped_WGS_reads/${sample}.sam --fast-local -p ${threads}"  >> ${pathScripts}/bsub_script_bowtie2
	 echo "${samtools} sort -@ ${threads} -o ${pathResults}/remapped_WGS_reads/${sample}.bam ${pathResults}/remapped_WGS_reads/${sample}.sam " >> ${pathScripts}/bsub_script_bowtie2
	 echo "bedtools genomecov -ibam ${pathResults}/remapped_WGS_reads/${sample}.bam -bg -split > ${pathResults}/remapped_WGS_reads/${sample}.bedGraph" >> ${pathScripts}/bsub_script_bowtie2

	 if [ ${cluster} = "pbil" ]; then
	     sbatch ${pathScripts}/bsub_script_bowtie2
	 fi
	 
	 if [ ${cluster} = "cloud" ]; then
	     chmod a+x ${pathScripts}/bsub_script_bowtie2
	     ${pathScripts}/bsub_script_bowtie2
	 fi
    done
fi

#########################################################################

if [ ${sp} = "Pauxi_pauxi" ]; then
    for sample in B10_S2 B10_S4 B51_S1 B51_S3
    do
	export pathR1=${pathWGS}/${sample}_R1_001_trimmed.fastq.gz
	export pathR2=${pathWGS}/${sample}_R2_001_trimmed.fastq.gz
	
	echo "#!/bin/bash" > ${pathScripts}/bsub_script_bowtie2
	
	 if [ ${cluster} = "pbil" ]; then
	     echo "#SBATCH --job-name=bowtie2_${sample}" >  ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --output=${pathScripts}/std_output_bowtie2_${sample}.txt" >>  ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --error=${pathScripts}/std_error_bowtie2_${sample}.txt" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --mem=5G" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_bowtie2
	     echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_bowtie2
	 fi

	 echo "bowtie2 -x ${pathResults}/${prefix} -1 ${pathR1} -2 ${pathR2} -S ${pathResults}/remapped_WGS_reads/${sample}.sam --fast-local -p ${threads}"  >> ${pathScripts}/bsub_script_bowtie2
	 echo "${samtools} sort -@ ${threads} -o ${pathResults}/remapped_WGS_reads/${sample}.bam ${pathResults}/remapped_WGS_reads/${sample}.sam " >> ${pathScripts}/bsub_script_bowtie2
	 echo "bedtools genomecov -ibam ${pathResults}/remapped_WGS_reads/${sample}.bam -bg -split > ${pathResults}/remapped_WGS_reads/${sample}.bedGraph" >> ${pathScripts}/bsub_script_bowtie2

	 if [ ${cluster} = "pbil" ]; then
	     sbatch ${pathScripts}/bsub_script_bowtie2
	 fi
	 
	 if [ ${cluster} = "cloud" ]; then
	     chmod a+x ${pathScripts}/bsub_script_bowtie2
	     ${pathScripts}/bsub_script_bowtie2
	 fi
    done
fi

#########################################################################
