#!/bin/bash

export library=$1
export cluster=$2
export threads=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
    export pathTools=/beegfs/home/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
    export pathTools=/sps/biometr/necsulea/Tools
fi

export pathWGS=${path}/data/WGS
export pathIndex=${path}/data/genome_indexes/MEGAHIT_RAGOUT
export pathResults=${path}/results/copy_number_variation
export pathScripts=${path}/scripts/copy_number_variation

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_process

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=${library}" >>  ${pathScripts}/bsub_script_process
    echo "#SBATCH --output=${pathScripts}/std_output_process_${library}.txt" >>  ${pathScripts}/bsub_script_process
    echo "#SBATCH --error=${pathScripts}/std_error_process_${library}.txt" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_process
fi

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=${library}" >>  ${pathScripts}/bsub_script_process
    echo "#SBATCH --output=${pathScripts}/std_output_process_${library}.txt" >>  ${pathScripts}/bsub_script_process
    echo "#SBATCH --error=${pathScripts}/std_error_process_${library}.txt" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_process
    echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_process
fi

#########################################################################

export run=0

#########################################################################
# Align the data

if [ -e ${pathResults}/${library}.bam ]; then
    echo "bwa already there"
else
    echo "bwa mem -t ${threads} -R \"@RG\tID:id\tSM:sample\tLB:lib\" ${pathIndex}/genome_sequence_renamed.fa ${pathWGS}/${library}_R1_001_trimmed.fastq.gz ${pathWGS}/${library}_R2_001_trimmed.fastq.gz \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > ${pathResults}/${library}.bam" >> ${pathScripts}/bsub_script_process

    export run=1
fi

#########################################################################
# Extract the discordant paired-end alignments.

if [ -e ${pathResults}/${library}.discordants.unsorted.bam ]; then
    echo "discordants already there"
else
    echo "samtools view -b -F 1294 ${pathResults}/${library}.bam > ${pathResults}/${library}.discordants.unsorted.bam" >> ${pathScripts}/bsub_script_process

    export run=1
fi

#########################################################################
# Extract the split-read alignments

if [ -e  ${pathResults}/${library}.splitters.unsorted.bam ]; then
    echo "split reads already there"
else
    echo "samtools view -h ${pathResults}/${library}.bam \
    | 	 ${pathTools}/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${pathResults}/${library}.splitters.unsorted.bam" >> ${pathScripts}/bsub_script_process

    export run=1
fi

#########################################################################
# Sort both alignments

if [ -e ${pathResults}/${library}.discordants ]; then
    echo "sorted discordants already there"
else
    echo "samtools sort ${pathResults}/${library}.discordants.unsorted.bam ${pathResults}/${library}.discordants" >> ${pathScripts}/bsub_script_process

    export run=1
fi

if [ -e ${pathResults}/${library}.splitters ]; then
    echo "sorted split reads already there"
else
    echo "samtools sort ${pathResults}/${library}.splitters.unsorted.bam ${pathResults}/${library}.splitters" >> ${pathScripts}/bsub_script_process

    export run=1
fi

#########################################################################
## Run everything

if [ ${run} = 1 ]; then
    if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
	sbatch ${pathScripts}/bsub_script_process
    fi

    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_process
	${pathScripts}/bsub_script_process
    fi
fi

#########################################################################
