#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3
export nthreads=$4

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

echo "#!/bin/bash" > ${pathScripts}/bsub_script_merge

#############################################################################

echo -n "samtools merge -o ${pathResults}/accepted_hits_all_samples.bam " >> ${pathScripts}/bsub_script_merge

for sample in `ls ${pathResults}`
do
    if [ -e ${pathResults}/${sample}/accepted_hits.bam ]; then
	echo -n ${pathResults}/${sample}/accepted_hits.bam" " >> ${pathScripts}/bsub_script_merge
    fi
done

echo "" >> ${pathScripts}/bsub_script_merge

#############################################################################

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_merge
    ${pathScripts}/bsub_script_merge
fi

#############################################################################
