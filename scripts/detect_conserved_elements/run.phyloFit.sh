#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export refsp=$1
export dataset=$2
export cluster=$3

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathAln=${path}/results/whole_genome_alignments/${dataset}/${refsp}
export pathSites=${path}/results/conserved_elements/${dataset}/${refsp}/4d/
export pathResults=${path}/results/conserved_elements/${dataset}/${refsp}/mod/

######################################################################

if [ -e ${pathResults} ]; then
    echo "output dir already there"
else
    mkdir -p ${pathResults}
fi

######################################################################

for chr in {1..29} Z W 
do
    export alnprefix="aln"

    if [ -e ${pathSites}/${alnprefix}_${chr}_4d-sites.ss ]; then
	
	if [ -e ${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod ]; then
	    echo "already done"
	else
	    phyloFit --tree ${pathAln}/tree.txt --msa-format SS --out-root ${pathResults}/phyloFit_nonconserved_4d-sites_${chr} ${pathSites}/${alnprefix}_${chr}_4d-sites.ss
	fi
    fi
done

######################################################################
