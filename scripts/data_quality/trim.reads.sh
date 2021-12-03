#!/bin/bash

export sp=$1
export sample=$2
export cluster=$3

###############################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi


if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

export pathWGS=${path}/data/WGS/${sp}

###############################################################

## This is cutadapt 2.8 with Python 3.8.10

if [ ${sp} = "Pauxi_pauxi" ]; then
    ## Nextera
    export adapter1=CTGTCTCTTATACACATCT
    export adapter2=CTGTCTCTTATACACATCT

    export pathR1=${pathWGS}/${sample}_R1_001.fastq.gz
    export pathR2=${pathWGS}/${sample}_R2_001.fastq.gz

    export pathR1out=${pathWGS}/${sample}_R1_001_trimmed.fastq.gz
    export pathR2out=${pathWGS}/${sample}_R2_001_trimmed.fastq.gz
fi

if [ ${sp} = "Basiliscus_vittatus" ]; then

    export pathR1=${pathWGS}/${sample}_R1.fastq.gz
    export pathR2=${pathWGS}/${sample}_R2.fastq.gz

    export pathR1out=${pathWGS}/${sample}_R1_trimmed.fastq.gz
    export pathR2out=${pathWGS}/${sample}_R2_trimmed.fastq.gz

    ## TruSeq
    export adapter1=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
    export adapter2=CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
fi

###############################################################

echo ${sample}

cutadapt --minimum-length 50 --trim-n -a ${adapter1} -A ${adapter2} -o ${pathR1out} -p ${pathR2out} ${pathR1} ${pathR2}

###############################################################
