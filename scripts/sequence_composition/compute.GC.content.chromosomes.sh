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

export pathEnsemblGenomes=${path}/data/genome_sequences/Ensembl103
export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
export pathScripts=${path}/scripts/sequence_composition

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export pathResults=${pathGenomeAssembly}
    export outfile=GCContent_ChromosomeSize.txt
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed.fa
    export pathResults=${pathGenomeAssembly}
    export outfile=GCContent_ChromosomeSize.txt
fi

#########################################################################

if [ ${assembly} = "Ensembl" ]; then
    export file=`ls ${pathEnsemblGenomes} | grep ${sp} | grep dna_sm.toplevel.fa.gz`
    export pathAssembly=${pathEnsemblGenomes}/${file}
    export pathResults=${pathEnsemblGenomes}
    export outfile=GCContent_ChromosomeSize_${sp}.txt
fi

#########################################################################

perl ${pathScripts}/compute.GC.content.chromosomes.pl --pathGenomeSequence=${pathAssembly} --pathOutput=${pathResults}/${outfile}

#########################################################################
