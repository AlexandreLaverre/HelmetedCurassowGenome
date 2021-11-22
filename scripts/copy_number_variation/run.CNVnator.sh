#!/bin/bash

set -e

export chr=$1
export cluster=$2

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathAln=${path}/results/bowtie2_alignments
export pathAssembly=${path}/results/genome_assembly/MEGAHIT_RAGOUT
export pathResults=${path}/results/CNVnator

##############################################################

## Extract read mapping

cnvnator -root ${pathResults}/allsamples_${chr}.root -tree ${pathAln}/accepted_hits_allsamples.bam -chrom ${chr}

##############################################################

# Generate histogram

cnvnator -root ${pathResults}/allsamples_${chr}.root -his 1000 -fasta ${pathAssembly}/genome_sequence_renamed.fa

##############################################################

# Calculate statistics

cnvnator -root ${pathResults}/allsamples_${chr}.root -stat 1000

##############################################################

# Partition

cnvnator -root ${pathResults}/allsamples_${chr}.root -partition 1000

##############################################################

# Call CNVs

cnvnator -root ${pathResults}/allsamples_${chr}.root -call 1000

##############################################################
