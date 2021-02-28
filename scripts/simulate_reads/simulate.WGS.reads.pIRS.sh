#!/bin/bash

export sp=$1
export prefixchr=$2

#########################################################################

export path=/beegfs/data/necsulea/IPLOSS
export pathGenome=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/simulated_WGS_reads/${sp}/pIRS
export pathScripts=${path}/scripts/simulate_reads

#########################################################################

## maximum read length 100bp to use the error rates already computed

if [ -e ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_100_500_1.fq.gz ]&&[ -e ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_100_500_2.fq.gz ]; then
    echo "already done"
else
    /beegfs/home/necsulea/Tools/pIRS_111/src/pirs/pirs simulate -i ${pathGenome}/${prefixchr}.fa -l 100 -x 30 -m 500 -v 25 -a 1 -g 1 -Q 33 -c 1 -o ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp
fi

#########################################################################

fastqc -o ${pathResults}/ ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_100_500_1.fq.gz ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_100_500_2.fq.gz

#########################################################################
