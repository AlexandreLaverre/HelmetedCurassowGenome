#!/bin/bash

export sp=$1
export trimmed=$2
export cluster=$3

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathData=${path}/data/RNASeq/${sp}

#############################################################################

## first combine all samples, paired_end only

if [ -e ${pathData}/allsamples_R1.fastq.gz ]||[ -e ${pathData}/allsamples_R2.fastq.gz ]; then
    echo "combined path exists"
else
    echo "combining samples"
    if [ ${trimmed} = "true" ]; then
	cat ${pathData}/*_R1_trimmed.fastq.gz > ${pathData}/allsamples_R1.fastq.gz
	cat ${pathData}/*_R2_trimmed.fastq.gz > ${pathData}/allsamples_R2.fastq.gz
    else
	cat ${pathData}/*_R1.fastq.gz > ${pathData}/allsamples_R1.fastq.gz
	cat ${pathData}/*_R2.fastq.gz > ${pathData}/allsamples_R2.fastq.gz
    fi
fi

#############################################################################
