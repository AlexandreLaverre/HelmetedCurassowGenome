#!/bin/bash

export prefix=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathHALParts=${path}/results/genome_assembly/MEGAHIT_RAGOUT/mafs_by_chr
export pathHAL=${path}/results/genome_assembly/MEGAHIT_RAGOUT

#########################################################################

export firstfile=`ls ${pathHALParts} | head -n 1`
echo ${firstfile}

cp ${pathHALParts}/${firstfile} ${pathHAL}/${prefix}.maf

for file in `ls ${pathHALParts} | grep maf | grep -v ${firstfile} `
do
    echo ${file}
    grep -v "#" ${pathHALParts}/${file} | sed '1d' >> ${pathHAL}/${prefix}.maf
done

#########################################################################
