#!/bin/bash

export library=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathWGS=${path}/data/WGS
export pathIndex=${path}/data/genome_indexes/MEGAHIT_RAGOUT
export pathResults=${path}/results/copy_number_variation
export pathScripts=${path}/scripts/copy_number_variation

#########################################################################
# Align the data

bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" ${pathIndex}/genome_sequence_renamed.fa ${pathWGS}/${library}_R1_001_trimmed.fastq.gz ${pathWGS}/${library}_R2_001_trimmed.fastq.gz \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > ${pathResults}/${library}.bam

#########################################################################
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 ${pathResults}/${library}.bam > ${pathResults}/${library}.discordants.unsorted.bam

#########################################################################

# Extract the split-read alignments
samtools view -h ${pathResults}/${library}.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${pathResults}/${library}.splitters.unsorted.bam

#########################################################################
# Sort both alignments
samtools sort ${pathResults}/${library}.discordants.unsorted.bam ${pathResults}/${library}.discordants
samtools sort ${pathResults}/${library}.splitters.unsorted.bam ${pathResults}/${library}.splitters

#########################################################################
