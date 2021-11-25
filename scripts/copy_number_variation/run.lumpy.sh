#!/bin/bash

export library=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

export pathResults=${path}/results/copy_number_variation
export pathScripts=${path}/scripts/copy_number_variation

#########################################################################

if [ ${library} = "all" ]; then
    export pathB=""
    export pathS=""
    export pathD=""

    for sample in B10_S2 B10_S4 B51_S3 B51_S1
    do
	export pathB=${pathResults}/${sample}.bam,${pathB}
	export pathS=${pathResults}/${sample}.splitters,${pathS}
	export pathD=${pathResults}/${sample}.discordants,${pathD}
    done
else
    export pathB=${pathResults}/${library}.bam
    export pathS=${pathResults}/${library}.splitters
    export pathD=${pathResults}/${library}.discordants
fi

#########################################################################

lumpyexpress -B ${pathB} -S ${pathS} -D ${pathD} -o ${pathResults}/${library}.vcf

#########################################################################
