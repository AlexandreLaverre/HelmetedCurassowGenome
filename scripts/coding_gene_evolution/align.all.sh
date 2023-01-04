#!/bin/bash

export nmax=$1
export step=$2
export cluster=$3
export dataset=$4

##########################################################################

export start=1


while [ ${start} -lt ${nmax} ]
do
    export end=$[$start+$step]

    echo ${start} ${end}
    
    ./align.coding.sequences.sh ${start} ${end} ${cluster} ${dataset}
    
    export start=$[$start+$step]
done

##########################################################################
