#!/bin/bash

export sp=$1
export prefixchr=$2
export method=$3

export readlen=150
export insertsize=500
export cov=30

#########################################################################

export path=/beegfs/data/necsulea/IPLOSS
export pathGenome=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/simulated_WGS_reads/${sp}/${method}
export pathScripts=${path}/scripts/simulate_reads

#########################################################################

bowtie2 -x ${pathGenome}/${prefixchr} -1 ${pathResults}/simulated_WGS_${prefixchr}_${cov}x_insertsize${insertsize}bp_${readlen}_${insertsize}_1.fq.gz -2 ${pathResults}/simulated_WGS_${prefixchr}_${cov}x_insertsize${insertsize}bp_${readlen}_${insertsize}_2.fq.gz -S ${pathResults}/simulated_WGS_${prefixchr}_${cov}x_insertsize${insertsize}bp_${readlen}_${insertsize}.${prefixchr}.sam --sensitive-local >& ${pathResults}/bowtie_output

#########################################################################
