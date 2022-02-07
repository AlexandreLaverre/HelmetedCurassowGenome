#!/bin/bash

export sp=$1
export cluster=$2
export nthreads=$3
export maxmem=$4

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathData=${path}/data/RNASeq/${sp}
export pathResults=${path}/results/transcriptome_assembly/${sp}/Trinity

## jellyfish 2.3.0
## samtools 1.10, htslib 1.10.2-3
## bowtie2 2.3.5.1
## salmon 1.6.0
## java openjdk 11.0.13 2021-10-19
## Trinity 2.13.2

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
    export libtype="--SS_lib_type RF"
else
    if [ ${sp} = "Chamaeleo_chamaeleon_recticrista" ]; then
	export libtype="" ## not stranded
    else
	echo "unknown library type!"
	exit
    fi
fi

#############################################################################

Trinity --seqType fq --max_memory ${maxmem} --left ${pathData}/allsamples_R1.fastq.gz --right ${pathData}/allsamples_R2.fastq.gz ${libtype} --CPU ${nthreads} --min_contig_length 200 --trimmomatic --output ${pathResults}

#############################################################################
