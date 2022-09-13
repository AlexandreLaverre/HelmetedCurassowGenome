#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export source=$3
export cluster=$4
export threads=$5
export hours=$6

#########################################################################

if [ ${cluster} = "pbil" ]; then
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

export speciesList=""
export pathsFastaProteins=""
export pathsTBlastNResults=""

for file in `ls ${pathResults} `
do
    export sp=`basename ${file} _vs_${assembly}.tblastn.out`
    export protfile=`ls ${pathProteinSequences} | grep "${sp}\." | grep fa`

    export speciesList=${sp},${speciesList}
    export pathsFastaProteins=${pathProteinSequences}/${protfile},${pathsFastaProteins}
    export pathsTBlastNResults=${pathResults}/${file},${pathsTBlastNResults}

    echo ${sp} ${protfile}
done   

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=orf_${target}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_${target}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_${target}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=orf_${refsp}" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --output=${pathScripts}/std_output_orf_${refsp}.txt" >>  ${pathScripts}/bsub_script_orf
    echo "#SBATCH --error=${pathScripts}/std_error_orf_${refsp}.txt" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orf
    echo "#SBATCH --time=${hours}:00:00" >> ${pathScripts}/bsub_script_orf
fi

echo "perl ${pathScripts}/extract.ORFs.pl --speciesList=${speciesList} --pathsFastaProteins=${pathsFastaProteins} --pathsTBlastNResults=${pathsTBlastNResults} --minORFLength=300 --minProteinFraction=0.25 --maxEValue=0.001 --maxGapFraction=0.1 --pathOutput=${pathResults}/CombinedORFs.txt">> ${pathScripts}/bsub_script_orf

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
    sbatch ${pathScripts}/bsub_script_orf
fi

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_orf
    ${pathScripts}/bsub_script_orf
fi

#########################################################################
