#!/bin/bash

set -e

export spset=$1
export dataset=$2
export phenannot=$3
export cluster=$4
export threads=$5
export multifilter=$6

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
    export pathTools=/beegfs/home/necsulea/Tools
fi

export pathResults=${path}/results/coding_gene_evolution/${spset}/${dataset}

##########################################################################

export OMP_NUM_THREADS=1

##########################################################################

## annotate tree

if [ -e ${pathResults}/pelican_annotated_tree_${phenannot}.txt ]; then
    echo "annotated tree already there"
else
    singularity exec -B ${path} ${pathTools}/pelican.sif pelican tree-annotation parsimony --output ${pathResults}/pelican_annotated_tree_${phenannot}.txt --tree ${pathResults}/species_tree_nobootstrap.txt --annotation ${pathResults}/phenotype_data_${phenannot}.txt
fi

##########################################################################

if [ "${multifilter}" != "1" ]; then
    singularity exec -B ${path} ${pathTools}/pelican.sif pelican scan discrete --multinomial-filter ${multifilter} --threads ${threads} --alphabet=nuc --name=${phenannot} --tree ${pathResults}/pelican_annotated_tree_${phenannot}.txt --alignment ${pathResults}/data_for_pelican --output ${pathResults}/pelican_output_${phenannot}_multinomial_${multifilter} --progress-bar
else
    singularity exec -B ${path} ${pathTools}/pelican.sif pelican scan discrete --threads ${threads} --alphabet=nuc --name=${phenannot} --tree ${pathResults}/pelican_annotated_tree_${phenannot}.txt --alignment ${pathResults}/data_for_pelican --output ${pathResults}/pelican_output_${phenannot} --progress-bar
fi

##########################################################################


