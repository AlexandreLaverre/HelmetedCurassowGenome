#!/bin/bash

export sp=$1
export cluster=$2
export nthreads=$3

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomes=${path}/data/genome_sequences/Ensembl103
export pathAssembly=${path}/results/genome_assembly/${sp}/MEGAHIT
export pathResults=${path}/results/genome_assembly/${sp}/MEGAHIT_RAGOUT

#########################################################################

export TMPDIR=${pathResults}/tmp

if [ -e ${TMPDIR} ]; then
    echo "tmpdir already there"
else
    echo "making directory"
    mkdir ${TMPDIR}
fi

#########################################################################

if [ ${sp} = "Pauxi_pauxi"]; then
    echo "(Anas_platyrhynchos_platyrhynchos, ((Penelope_pileata, Pauxi_pauxi), Gallus_gallus));" > ${pathResults}/seqFile
    echo "" >>  ${pathResults}/seqFile
    
    for sp in Anas_platyrhynchos_platyrhynchos Penelope_pileata  Gallus_gallus
    do
	echo "${sp} ${pathGenomes}/${sp}/genome_sm_clean.fa" >> ${pathResults}/seqFile
    done
    
    echo "Pauxi_pauxi ${pathAssembly}/final.contigs.clean.fa"  >> ${pathResults}/seqFile
fi

#########################################################################

if [ ${sp} = "Basiliscus_vittatus"]; then
    echo "(Salvator_merianae, ((Basiliscus_vittatus, Anolis_carolinensis), Pseudonaja_textilis));" > ${pathResults}/seqFile
    echo "" >>  ${pathResults}/seqFile
    
    for sp in Salvator_merianae Anolis_carolinensis Pseudonaja_textilis
    do
	echo "${sp} ${pathGenomes}/${sp}/genome_sm_clean.fa" >> ${pathResults}/seqFile
    done

    if [ -e ${pathAssembly}/final.contigs.clean.fa ]; then
	echo "clean file already there"
    else
	perl ${pathScripts}/cleanup.fasta.names.pl --pathInput=${pathAssembly}/final.contigs.fa --pathOutput=${pathAssembly}/final.contigs.clean.fa
    fi
    
    echo "Basiliscus_vittatus ${pathAssembly}/final.contigs.clean.fa"  >> ${pathResults}/seqFile
fi

#########################################################################

docker run -v ${path}:/mnt/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 cactus --binariesMode local --workDir ${pathResults}/ --maxServiceJobs 1 --maxCores ${nthreads} --maxMemory 60G ${pathResults}/jobStore ${pathResults}/seqFile ${pathResults}/alignment.hal

#########################################################################
