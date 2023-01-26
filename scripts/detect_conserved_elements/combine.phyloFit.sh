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
    export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved} 
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_macro_chromosomes.mod

######################################################################

## micro-chromosomes

export pathsNonConserved=""

for chr in {9..16} {18..29} 
do
    export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved} 
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_micro_chromosomes.mod

######################################################################

## sex-chromosomes (simple cp for phast pipeline homogenisation)

cp ${pathResults}/phyloFit_nonconserved_4d-sites_Z.mod ${pathResults}/phyloFit_nonconserved_4d-sites_avg_sex_chromosomes.mod

######################################################################

## micro/macro/sex (for scaffold)

export pathsNonConserved=""

for chr in {1..16} {18..29} Z
do
    export pathsNonConserved="${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod",${pathsNonConserved} 
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_all_chromosomes.mod

######################################################################


