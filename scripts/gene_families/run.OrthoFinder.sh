#!/bin/bash

export type=$1  ## tree inference method
export cluster=$2
export threads=$3

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

export pathResults=${path}/results/gene_families/OrthoFinder
export pathScripts=${path}/scripts/gene_families

##########################################################################

## OrthoFinder 2.5.4
## MAFFT v7.453
## mmqseqs 2 
## IQ-TREE multicore version 1.6.12 for Linux 64-bit 

##########################################################################

if [ ${cluster} = "cloud" ]; then
    ulimit -n 50000
fi

##########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_orthofinder

##########################################################################

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=orthofinder_${type}" >>  ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --output=${pathScripts}/std_output_orthofinder.txt" >>  ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --error=${pathScripts}/std_error_orthofinder.txt" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --partition=bigmem" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --mem=64G" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --time=168:00:00" >> ${pathScripts}/bsub_script_orthofinder
fi
    
##########################################################################

echo "orthofinder -f ${pathResults} -o ${pathResults}/${type} -t ${threads} -a ${threads} -I 2 -S mmseqs -A muscle -M msa -y -T ${type}" >>  ${pathScripts}/bsub_script_orthofinder

##########################################################################

if [ ${cluster} = "pbil" ]; then
    sbatch ${pathScripts}/bsub_script_orthofinder
fi

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_orthofinder
    ${pathScripts}/bsub_script_orthofinder
fi

##########################################################################
