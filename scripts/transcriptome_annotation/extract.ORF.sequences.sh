#!/bin/bash

########################################################################

export target=$1
export cluster=$2

export assembly="Trinity"

#########################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbildeb" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

export pathResults=${path}/results/transcriptome_assembly/${target}
export pathScripts=${path}/scripts/transcriptome_annotation

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=orf_${ref}_${target}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_seq.txt_${target}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_seq_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=4G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=1:00:00" >> ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=orf_${refsp}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_seq_${target}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_seq_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=4G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=1:00:00" >> ${pathScripts}/bsub_script_orf
fi

echo "perl ${pathScripts}/extract.ORF.sequences.pl  --pathORFs=${pathResults}/CombinedORFs.txt --pathContigs=${pathResults}/Trinity.fasta --pathGeneticCode=${pathScripts}/standard_genetic_code.txt --pathOutputCDS=${pathResults}/CombinedORFs_CDS.fa  --pathOutputProtein=${pathResults}/CombinedORFs_Proteins.fa">> ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
    sbatch ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbildeb" ]; then
    chmod a+x ${pathScripts}/bsub_script_orf
    ${pathScripts}/bsub_script_orf
fi

#########################################################################
