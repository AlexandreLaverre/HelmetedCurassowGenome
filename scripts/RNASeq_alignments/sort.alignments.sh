#!/bin/bash

export sp=$1
export assembly=$2
export sample=$3
export cluster=$4
export nthreads=$5

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/RNASeq_alignments/${sp}/${assembly}
export pathScripts=${path}/scripts/RNASeq_alignments

#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.bam ]; then
    echo "already done"
    exit
fi

#############################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_sort

#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.sam ]; then
    echo "samtools sort -m 5G -@ ${nthreads} -o ${pathResults}/${sample}/accepted_hits.bam -O bam ${pathResults}/${sample}/accepted_hits.sam" >> ${pathScripts}/bsub_script_sort
else
    echo "cannot find sam file"
fi

#############################################################################

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_sort
    ${pathScripts}/bsub_script_sort
fi

#############################################################################
