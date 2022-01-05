#!/bin/bash

export sp=$1
export cluster=$2

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathData=${path}/data/RNASeq/${sp}
export pathResults=${path}/results/transcriptome_assembly/${sp}/TransLiG

## TransLiG v1.3 

#############################################################################

if [ -e ${pathResults} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}
fi

#############################################################################

## first combine all samples, paired_end only

if [ -e ${pathData}/allsamples_R1.fastq.gz ]||[ -e ${pathData}/allsamples_R2.fastq.gz ]; then
    echo "combined path exists"
else
    echo "combining samples"
    cat ${pathData}/*_R1_trimmed.fastq.gz > ${pathData}/allsamples_R1.fastq.gz
    cat ${pathData}/*_R2_trimmed.fastq.gz > ${pathData}/allsamples_R2.fastq.gz
fi

#############################################################################

## RF: TruSeq mRNA stranded, also Kapa Biosystems Stranded mRNA

if [ ${sp} = "Chamaeleo_calyptratus" ]; then
    export libtype="RF"
else
    echo "unknown library type!"
    exit
fi

#############################################################################

TransLiG -m ${libtype} -s fq -p pair -l ${pathData}/allsamples_R1.fastq.gz -r ${pathData}/allsamples_R2.fastq.gz -o ${pathResults}

#############################################################################

