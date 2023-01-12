#!/bin/bash

set -e

export spset=$1
export dataset=$2
export phenannot=$3
export cluster=$4
export threads=$5

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

singularity exec -B ${path} ${pathTools}/pelican.sif pelican tree-parsimony-annotation --output ${pathResults}/pelican_annotated_tree_${phenannot}.txt --tree ${pathResults}/species_tree_nobootstrap.txt --annotation ${pathResults}/phenotype_data_${phenannot}.txt

##########################################################################

singularity exec -B ${path} ${pathTools}/pelican.sif pelican scan --multinomial-filter 0.001 --threads ${threads} --translate --tree ${pathResults}/pelican_annotated_tree_${phenannot}.txt --alignment ${pathResults}/data_for_pelican --output ${pathResults}/pelican_output_${phenannot} --progress-bar

##########################################################################


