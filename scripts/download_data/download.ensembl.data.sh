#!/bin/bash

export cluster=$1

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathProteins=${path}/data/protein_sequences/Ensembl103

##########################################################################

for file in `ls ${pathProteins} | grep pep.all.fa`
do
    export prefix=`basename ${file} .pep.all.fa`
    export sp=`echo ${prefix} | cut -f 1 -d '.'`
    export lowersp=`echo "$sp" | awk '{print tolower($0)}'`

    echo ${file} ${prefix} ${sp} ${lowersp}

    wget http://ftp.ensembl.org/pub/release-103/gff3/${lowersp}/${sp}.${prefix}.gff3.gz
done

##########################################################################

