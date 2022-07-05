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

if [ ${cluster} = "pbil_local" ]; then
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

echo "#!/bin/bash" > ${pathScripts}/bsub_script_bowtie2

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=bowtie2_index" >>  ${pathScripts}/bsub_script_bowtie2
    echo "#SBATCH --output=${pathScripts}/std_output_bowtie2_index.txt" >>  ${pathScripts}/bsub_script_bowtie2
    echo "#SBATCH --error=${pathScripts}/std_error_bowtie2_index.txt" >> ${pathScripts}/bsub_script_bowtie2
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_bowtie2
    echo "#SBATCH --mem=5G" >> ${pathScripts}/bsub_script_bowtie2
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_bowtie2
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_bowtie2
fi

echo "bowtie2-build --threads ${threads} ${pathAssembly} ${pathResults}/${prefix}">> ${pathScripts}/bsub_script_bowtie2

if [ ${cluster} = "pbil" ]; then
    sbatch ${pathScripts}/bsub_script_bowtie2
fi

if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbil_local" ]; then
    chmod a+x ${pathScripts}/bsub_script_bowtie2
    ${pathScripts}/bsub_script_bowtie2
fi


#########################################################################
