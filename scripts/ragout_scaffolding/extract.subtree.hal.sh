#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export pathIPLOSS=/mnt/mydatalocal/IPLOSS
    export pathHCG=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${pathIPLOSS}/data/whole_genome_alignments/363birds
export pathResults=${pathHCG}/data/whole_genome_alignments/363birds

#########################################################################

halExtract --root birdAnc341 ${pathHAL}/363-avian-2020.hal ${pathResults}/galloanserae.hal 

#########################################################################
