#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export refsp=$3
export source=$4
export cluster=$5
export threads=$6
export hours=$7

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

export pathProteinSequences=${path}/data/protein_sequences/${source}/primary_transcripts
export pathTranscriptomeAssembly=${path}/results/transcriptome_assembly/${target}
export pathResults=${path}/results/transcriptome_assembly/${target}/tblastn_results
export pathScripts=${path}/scripts/transcriptome_annotation

#########################################################################

export protfile=`ls ${pathProteinSequences} | grep "${refsp}\." | grep fa`

export pathFastaProteins=${pathProteinSequences}/${protfile}
export pathTBlastNResults=${pathResults}/${refsp}_vs_${assembly}.tblastn.out

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=orf_${ref}_${target}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_${ref}_${target}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_${ref}_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=orf_${refsp}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_${refsp}_${target}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_${refsp}_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_orf
fi

echo "perl ${pathScripts}/extract.ORFs.pl  --pathFastaProteins=${pathFastaProteins} --pathTBlastNResults=${pathTBlastNResults} --minORFLength=300 --minProteinFraction=0.25 --maxEValue=0.001 --maxGapFraction=0.1 --pathOutput=${pathResults}/${refsp}_ORFs.txt">> ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
    sbatch ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbildeb" ]; then
    chmod a+x ${pathScripts}/bsub_script_orf
    ${pathScripts}/bsub_script_orf
fi

#########################################################################
