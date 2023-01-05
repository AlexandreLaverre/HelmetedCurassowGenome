#!/bin/bash

########################################################################

export target=$1
export cluster=$2
export threads=$3
export hours=$4

#########################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbildeb" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi


export pathResults=${path}/results/transcriptome_assembly/${target}
export pathScripts=${path}/scripts/transcriptome_annotation

export prefix=CombinedORFs_Proteins

#########################################################################

perl ${pathScripts}/extract.ORF.clusters.pl --pathInputCDS=${pathResults}/CombinedORFs_CDS.fa --pathInputProteins=${pathResults}/CombinedORFs_Proteins.fa --pathDiamondResults=${pathResults}/CombinedORFs_Proteins_vs_self.diamond.blastp.out --minPercentageIdentity=95 --minAlignedLengthFraction=0.75 --pathOutputClusters=${pathResults}/CombinedORFs_Clusters.txt --pathOutputCDS=${pathResults}/CombinedORFs_RepresentativeORF_CDS.fa  --pathOutputProteins=${pathResults}/CombinedORFs_RepresentativeORF_Proteins.fa 

#########################################################################
