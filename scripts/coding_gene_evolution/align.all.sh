#!/bin/bash

export nmax=$1
export cluster=$2

##########################################################################

export start=1


while [ ${start} -lt ${nmax} ]
do
    export end=$[$start+1000]

    echo ${start} ${end}
    
    ./align.coding.sequences.sh ${start} ${end} ${cluster}
    
    export start=$[$start+1000]
done

##########################################################################
