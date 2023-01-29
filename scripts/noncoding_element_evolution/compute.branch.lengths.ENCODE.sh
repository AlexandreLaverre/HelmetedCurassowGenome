#!/bin/bash

set -e 

export cluster=$1
export method=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${path}/results/whole_genome_alignments/avian_366
export pathResults=${path}/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse
export pathScripts=${path}/scripts/noncoding_element_evolution

#########################################################################

if [ -e ${pathResults}/${method}_results ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}/${method}_results
fi

#########################################################################

if [ -e ${pathResults}/aln_by_element ]; then
    for file in `ls ${pathResults}/aln_by_element | grep filtered.phy$`
    do
	export prefix=`basename ${file} .filtered.phy`

	if [ -e ${pathResults}/${method}_results/${prefix}.start.tree ]; then
	    echo "tree already there"
	else
	    Rscript --vanilla extract.species.tree.R ${pathHAL}/Birds366_tree.txt ${pathResults}/aln_by_element/${prefix}.filtered.fa ${pathResults}/${method}_results/${prefix}.start.tree
	fi
       
	if [ ${method} = "iqtree" ]; then
	    if [ -e ${pathResults}/iqtree_results/${prefix}.treefile ]; then
		echo "already done"
	    else
		iqtree -s ${pathResults}/aln_by_element/${prefix}.filtered.fa -st DNA  -te ${pathResults}/iqtree_results/${prefix}.start.tree -m GTR --pre ${pathResults}/iqtree_results/${prefix} 
	    fi
	fi

	if [ ${method} = "phyml" ]; then
	    if [ -e ${pathResults}/phyml_results/${prefix}.filtered.phy_phyml_tree.txt ]; then
		echo "already done"
	    else
		phyml -i ${pathResults}/aln_by_element/${prefix}.filtered.phy -d nt --sequential -m GTR -o lr  -u ${pathResults}/phyml_results/${prefix}.start.tree --no_memory_check --quiet
		mv ${pathResults}/aln_by_element/${prefix}.filtered.phy*phyml*  ${pathResults}/phyml_results
	    fi
	fi

    done
fi

#########################################################################
