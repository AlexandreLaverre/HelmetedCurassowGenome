#!/bin/bash

export sp=$1
export prefixchr=$2

#########################################################################

export path=/beegfs/data/necsulea/IPLOSS
export pathGenome=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/simulated_WGS_reads/${sp}/ART
export pathScripts=${path}/scripts/simulate_reads

#########################################################################

if [ -e ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_1.fq.gz ]&&[ -e ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_2.fq.gz ]; then
    echo "already done"
else

    art_illumina -ss HS25 -i ${pathGenome}/${prefixchr}.fa -p -l 150 -f 30 -m 500 -s 25 -o ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_
    
    gzip ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_*fq
fi

#########################################################################

fastqc -o ${pathResults}/ ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_1.fq.gz ${pathResults}/simulated_WGS_${prefixchr}_30x_insertsize500bp_150_500_2.fq.gz

#########################################################################
