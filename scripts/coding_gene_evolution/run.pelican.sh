#!/bin/bash

set -e

export dataset=$1
export cluster=$2
export threads=$3

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools
fi

export pathResults=${path}/results/coding_gene_evolution/

##########################################################################

export OMP_NUM_THREADS=1

##########################################################################

## annotate tree

singularity exec -B ${path} ${pathTools}/pelican.sif pelican tree-parsimony-annotation --output ${pathResults}/pelican_annotated_tree_${dataset}.txt --tree ${pathResults}/species_tree_nobootstrap.txt --annotation ${pathResults}/phenotype_data_${dataset}.txt

##########################################################################

singularity exec -B ${path} ${pathTools}/pelican.sif pelican scan --multinomial-filter 0.001 --threads ${threads} --translate --tree ${pathResults}/pelican_annotated_tree_${dataset}.txt --alignment ${pathResults}/data_for_pelican --output ${pathResults}/pelican_output_${dataset} --progress-bar

##########################################################################


