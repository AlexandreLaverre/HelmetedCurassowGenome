#!/bin/bash

export start=$1
export end=$2
export cluster=$3

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathCDS=${path}/data/coding_sequences
export pathAnnot=${path}/results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined
export pathResults=${path}/results/coding_gene_evolution
export pathScripts=${path}/scripts/coding_gene_evolution

##########################################################################

echo "#!/bin/bash" > ${pathScripts}/log/bsub_prank_${start}_${end}

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=prank_${start}_${end}" >> ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --partition=normal" >>${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --output=${pathScripts}/log/std_out_prank_${start}_${end}" >> ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --error=${pathScripts}/log/std_err_prank_${start}_${end}" >> ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --time=2:00:00" >>  ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --mem=20G" >> ${pathScripts}/log/bsub_prank_${start}_${end}
fi


if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=prank_${start}_${end}" >> ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --output=${pathScripts}/log/std_out_prank_${start}_${end}" >> ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --error=${pathScripts}/log/std_err_prank_${start}_${end}" >> ${pathScripts}/log/bsub_prank_${start}_${end}
    echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/log/bsub_prank_${start}_${end}
fi

##########################################################################

for file in `ls ${pathResults}/CDS | grep unaln | sed -n "${start},${end}p;${end}q"`
do
    export prefix=`basename ${file} .unaln.fa`

    echo ${prefix}
    
    echo "prank -codon -once -f=fasta -t=${pathResults}/CDS/${prefix}.tree -d=${pathResults}/CDS/${file} -o=${pathResults}/CDS/${prefix}.aln.fa" >> ${pathScripts}/log/bsub_prank_${start}_${end}

    echo "prank -convert -f=phylips -d=${pathResults}/CDS/${prefix}.aln.fa -o=${pathResults}/CDS/${prefix}.aln.phy"  >> ${pathScripts}/log/bsub_prank_${start}_${end}
    
done

##########################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
    sbatch ${pathScripts}/log/bsub_script_prank_${start}_${end}
fi

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/log/bsub_script_prank_${start}_${end}
    ${pathScripts}/log/bsub_script_prank_${start}_${end}
fi

##########################################################################
