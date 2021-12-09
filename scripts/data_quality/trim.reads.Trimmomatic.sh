#!/bin/bash

export sp=$1
export datatype=$2
export seqtype=$3
export cluster=$4
export threads=$5

###############################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

export pathData=${path}/data/${datatype}/${sp}

## trimmomatic/focal,now 0.39+dfsg-1

###############################################################

if [ ${USER} = "alaverre" ]; then
    export pathAdapt=/beegfs/data/alaverre/Tools/envs/share/trimmomatic/adapters
    export trimmomatic="trimmomatic "
else
    export pathAdapt=/usr/share/trimmomatic
    export trimmomatic="Trimmomatic"
fi

###############################################################

if [ ${sp} = "Pauxi_pauxi" ]; then
    export sequence=NexteraPE-PE.fa
fi

if [ ${sp} = "Basiliscus_vittatus" ]; then
    export sequence=TruSeq2-PE.fa
fi

if [ ${sp} = "Chamaeleo_calyptratus" ]; then
    export sequence=TruSeq3-PE-2.fa
fi

###############################################################

if [ ${seqtype} = "single_end" ]; then
    for path in `ls ${pathData} | grep fastq.gz | grep -v trimmed | grep -v R1 | grep -v R2`
    do
	export file=`basename ${path} .fastq.gz`

	if [ -e ${pathData}/${file}_trimmed.fastq.gz ]; then
	    echo ${file}" already done"
	else
	    ${trimmomatic}SE -phred33 -threads ${threads} ${pathData}/${path} ${pathData}/${file}_trimmed.fastq.gz ILLUMINACLIP:/${pathAdapt}/${sequence}:2:40:15
	fi
    done
fi

###############################################################

if [ ${seqtype} = "paired_end" ]; then
   
    for path in `ls ${pathData} | grep fastq.gz | grep -v trimmed | grep R2 | grep fastq.gz `
    do
	export file=`basename ${path} _R2.fastq.gz`

	if [ -e ${pathData}/${file}_R1_trimmed.fastq.gz ]&&[ -e ${pathData}/${file}_R2_trimmed.fastq.gz ]; then
	    echo ${file}" already done"
	else
	    ${trimmomatic}PE -phred33 -threads ${threads} ${pathData}/${file}_R1.fastq.gz ${pathData}/${file}_R2.fastq.gz ${pathData}/${file}_R1_trimmed.fastq.gz ${pathData}/${file}_R1_trimmed_unpaired.fastq.gz ${pathData}/${file}_R2_trimmed.fastq.gz ${pathData}/${file}_R2_trimmed_unpaired.fastq.gz ILLUMINACLIP:/${pathAdapt}/${sequence}:2:40:15
	fi
    done
fi

###############################################################
