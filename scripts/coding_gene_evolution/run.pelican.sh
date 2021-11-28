#!/bin/bash

export cluster=$1
export threads=$2

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/home/ubuntu/data/mydatalocal/Tools
fi

export pathResults=${path}/results/coding_gene_evolution/

##########################################################################

export OMP_NUM_THREADS=1

##########################################################################

## annotate tree

singularity exec -B ${path} ${pathTools}/pelican.sif pelican tree-parsimony-annotation --output ${pathResults}/pelican_annotated_tree.txt --tree ${pathResults}/species_tree_nobootstrap.txt --annotation ${pathResults}/phenotype_data.txt

##########################################################################

singularity exec -B ${path} ${pathTools}/pelican.sif pelican scan ---multinomial-filter 0.001 --threads ${threads} --translate --tree ${pathResults}/pelican_annotated_tree.txt --alignment ${pathResults}/data_for_pelican --output ${pathResults}/pelican_output --progress-bar

##########################################################################


