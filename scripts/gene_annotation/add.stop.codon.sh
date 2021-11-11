#!/bin/bash

export ref=$1
export source=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathGenomes=${path}/data/genome_sequences/${source}
export pathAnnotations=${path}/data/genome_annotations/${source}
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

export genomefile=`ls ${pathGenomes} | grep ${ref}'\.' | grep fa | grep -v fai$`
export annotfile=`ls ${pathAnnotations} | grep ${ref}'\.' | grep gff`

#########################################################################

echo "genome file "${genomefile}
echo "annot file "${annotfile}

#########################################################################

perl ${pathScripts}/add.stop.codon.pl --pathInputGFF=${pathAnnotations}/${annotfile} --pathOutputGFF=${pathAnnotations}/${ref}_withstopcodons.gff

# -J   discard any mRNAs that either lack initial START codon
#   or the terminal STOP codon, or have an in-frame stop codon
#   (i.e. only print mRNAs with a complete CDS)

gffread -J -g ${pathGenomes}/${genomefile} ${pathAnnotations}/${ref}_withstopcodons.gff > ${pathAnnotations}/${ref}_filtered.gff

#########################################################################
