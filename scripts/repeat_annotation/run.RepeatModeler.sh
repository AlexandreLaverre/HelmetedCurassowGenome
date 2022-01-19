#!/bin/bash

########################################################################

export sp=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi


export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathResults=${path}/results/repeats/${assembly}/${sp}/RepeatModeler
export pathScripts=${path}/scripts/repeat_annotation

# Search Engine = rmblast 2.11.0+
# Dependencies: TRF 4.09, RECON , RepeatScout 1.0.6, RepeatMasker 4.1.2
# LTR Structural Analysis: Enabled ( GenomeTools 1.6.1, LTR_Retriever v2.9.0,
#                                    Ninja 0.95-cluster_only, MAFFT 7.475,
#                                    CD-HIT 4.8.1 )

#########################################################################

RepeatModeler -database ${pathResults}/repeat_modeler_db -pa 8 -LTRStruct >& ${pathResults}/RepeatModeler.out

#########################################################################
