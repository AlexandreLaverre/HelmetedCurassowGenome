#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomes=${path}/data/genome_sequences
export pathAssembly=${path}/results/genome_assembly/MEGAHIT
export pathResults=${path}/results/genome_assembly/MEGAHIT_RAGOUT

#########################################################################

export TMPDIR=${pathResults}/tmp

if [ -e ${TMPDIR} ]; then
    echo "tmpdir already there"
else
    echo "making directory"
    mkdir ${TMPDIR}
fi

#########################################################################

echo "(Anas_platyrhynchos_platyrhynchos, ((Penelope_pileata, Pauxi_pauxi), Gallus_gallus));" > ${pathResults}/seqFile
echo "" >>  ${pathResults}/seqFile

for sp in Anas_platyrhynchos_platyrhynchos Penelope_pileata  Gallus_gallus
do
    echo "${sp} ${pathGenomes}/${sp}/genome_sm_clean.fa" >> ${pathResults}/seqFile
done

echo "Pauxi_pauxi ${pathAssembly}/final.contigs.clean.fa"  >> ${pathResults}/seqFile

#########################################################################

docker run -v ${path}:/mnt/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 cactus --maxServiceJobs 4 --maxCores 12 --maxMemory 50G --maxDisk 600G --defaultDisk 20G --defaultMemory 5G --defaultCores 4  --binariesMode local --workDir ${pathResults}/ ${pathResults}/jobStore ${pathResults}/seqFile ${pathResults}/alignment.hal

#########################################################################
