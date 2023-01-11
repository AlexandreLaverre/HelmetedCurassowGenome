#!/bin/bash

export nmax=$1
export step=$2
export cluster=$3
export spset=$4
export dataset=$5

##########################################################################

export start=1


while [ ${start} -lt ${nmax} ]
do
    export end=$[$start+$step]

    echo ${start} ${end}
    
    ./align.coding.sequences.sh ${start} ${end} ${cluster} ${spset} ${dataset}
    
    export start=$[$start+$step]
done

##########################################################################
