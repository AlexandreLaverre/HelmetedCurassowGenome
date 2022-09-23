#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export cluster=$3
export threads=$4
export hours=$5

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

export pathTranscriptomeAssembly=${path}/results/transcriptome_assembly/${target}
export pathResults=${path}/results/transcriptome_assembly/${target}/tblastn_results
export pathScripts=${path}/scripts/transcriptome_annotation

#########################################################################

export pathsORFs=""

for file in `ls ${pathResults} | grep ORF | grep -v Combined`
do
    export pathsORFs=${pathResults}/${file},${pathsORFs}
done

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=orf_${ref}_${target}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_${ref}_${target}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_${ref}_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=orf_${refsp}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_${refsp}_${target}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_${refsp}_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_orf
fi

echo "perl ${pathScripts}/combine.ORFs.pl  --pathsORFs=${pathsORFs} --pathOutput=${pathResults}/CombinedORFs.txt">> ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
    sbatch ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbildeb" ]; then
    chmod a+x ${pathScripts}/bsub_script_orf
    ${pathScripts}/bsub_script_orf
fi

#########################################################################
