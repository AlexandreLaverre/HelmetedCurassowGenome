#!/bin/bash

export assembly=$1
export ref=$2
export source=$3
export cluster=$4
export cores=$5

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathResults=${path}/results/genome_annotation/${assembly}/genBlastG/${ref}
export pathProteins=${path}/data/protein_sequences/${source}

if [ ${source} = "Ensembl103" ]; then
    export pathProteins=${path}/data/protein_sequences/${source}/primary_transcripts
fi

#########################################################################

if [ -e ${pathResults} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
    export suffix=final.contigs
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
    export suffix=genome_sequence
fi

#########################################################################

export refprotfile=`ls ${pathProteins} | grep ${ref} | grep pep.all.fa`

#########################################################################

if [ ${cluster} = "cloud" ]; then
    cd /mnt/mydatalocal/Tools/genBlast_v138_linux_x86_64
fi

genblast -p genblastg -t ${pathAssembly} -q ${pathProteins}/${refprotfile} -e 1e-2 -g T -f F -a 0.5 -d 100000 -r 10 -c 0.5 -s 0 -i 15 -x 20 -n 20 -v 2 -h 0 -j 3 -norepair -gff -cdna -pro -o ${pathResults}/annotation

#########################################################################
