#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export refsp=$1
export dataset=$2
export cluster=$3

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/conserved_elements/${dataset}/${refsp}/mod/

######################################################################

## macro-chromosomes

export pathsNonConserved=""

for chr in {1..8}
do
    if [ -e ${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod ]; then
	export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved}
    fi
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_macro_chromosomes.mod

######################################################################

## micro-chromosomes

export pathsNonConserved=""

for chr in {9..33} 
do
    if [ -e ${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod ]; then
	export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved}
    fi
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_micro_chromosomes.mod

######################################################################

## sex-chromosomes 

export pathsNonConserved=""

for chr in Z W
do
    if [ -e ${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod ]; then
	export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved}
    fi
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_sex_chromosomes.mod


######################################################################

## micro/macro/sex (for scaffold)

export pathsNonConserved=""

for chr in {1..33} Z
do
    if [ -e ${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod ]; then
	export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved}
    fi
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_all_chromosomes.mod

######################################################################


