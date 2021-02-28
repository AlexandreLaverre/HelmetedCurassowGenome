#!/bin/bash

export sp=$1
export prefixchr=$2

#########################################################################

export path=/beegfs/data/necsulea/IPLOSS
export pathGenome=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/simulated_WGS_reads/${sp}/wgsim
export pathScripts=${path}/scripts/simulate_reads

#########################################################################

if [ -e ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_1.fq.gz ]&&[ -e ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_2.fq.gz ]; then
    echo "already done"
else
    wgsim -e 0.0005 -1 150 -2 150 -N 55000000 -d 500 -s 25 ${pathGenome}/${prefixchr}.fa ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_1.fq ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_2.fq
    
    gzip ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_1.fq
    gzip ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_2.fq
fi

#########################################################################

fastqc -o ${pathResults}/ ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_1.fq.gz ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_2.fq.gz

#########################################################################
