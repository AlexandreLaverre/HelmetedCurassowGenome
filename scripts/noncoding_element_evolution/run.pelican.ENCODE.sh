#!/bin/bash

set -e

export cluster=$1
export threads=$2
export phenannot=$3
export multifilter=$4

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
    export pathTools=/beegfs/home/necsulea/Tools
fi

export pathResults=${path}/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse
export pathResultsCoding=${path}/results/coding_gene_evolution/all_species/birds

##########################################################################

export OMP_NUM_THREADS=1

##########################################################################

cp ${pathResultsCoding}/pelican_annotated_tree_${phenannot}.txt ${pathResults}

##########################################################################

if [ -e ${pathResults}/data_for_pelican ]; then
    echo "input files already there"
else
    mkdir -p ${pathResults}/data_for_pelican
    
    for file in `ls ${pathResults}/aln_by_element/ | grep filtered.fa`
    do
	cp ${pathResults}/aln_by_element/${file} ${pathResults}/data_for_pelican/
    done
fi

##########################################################################

singularity exec -B ${path} ${pathTools}/pelican.sif pelican scan discrete --multinomial-filter ${multifilter} --threads ${threads} --alphabet=AA --name=${phenannot} --tree ${pathResults}/pelican_annotated_tree_${phenannot}.txt --alignment ${pathResults}/data_for_pelican --output ${pathResults}/pelican_output_${phenannot}_multinomial_${multifilter} --progress-bar

##########################################################################


