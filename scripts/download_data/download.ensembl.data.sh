#!/bin/bash

export cluster=$1

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathProteins=${path}/data/protein_sequences/Ensembl103
export pathAnnot=${path}/data/genome_annotations/Ensembl103
export pathGenomes=${path}/data/genome_sequences/Ensembl103

###########################################################################

for file in `ls ${pathProteins} | grep pep.all.fa`
do
    export prefix=`basename ${file} .pep.all.fa`
    export sp=`echo ${prefix} | cut -f 1 -d '.'`
    export lowersp=`echo "$sp" | awk '{print tolower($0)}'`

    echo ${file} ${prefix} ${sp} ${lowersp}
    
    ## gff3 files
    wget http://ftp.ensembl.org/pub/release-103/gff3/${lowersp}/${prefix}.103.gff3.gz

    gunzip ${prefix}.103.gff3.gz
    
    mv ${prefix}.103.gff3 ${pathAnnot}

    ## genome sequence
    wget http://ftp.ensembl.org/pub/release-103/fasta/${lowersp}/dna/${prefix}.dna.toplevel.fa.gz

    gunzip ${prefix}.dna.toplevel.fa.gz

    mv ${prefix}.dna.toplevel.fa ${pathGenomes}
done

##########################################################################



