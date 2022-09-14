#!/bin/bash

########################################################################

export target=$1
export assembly=$2
export source=$3
export cluster=$4

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
export pathGenomeAnnotation=${path}/results/genome_annotation/${target}/${assembly}/GeMoMa
export pathResults=${pathGenomeAnnotation}/combined/diamond_results
export pathScripts=${path}/scripts/filter_gene_annotation

#########################################################################

export speciesList=""
export pathsFastaProteins=""
export pathsDiamondResults=""

for sp in `ls ${pathGenomeAnnotation} | grep -v combined | grep -v "\." `
do
    export file=GeMoMa_vs_${sp}.diamond.blastp.out
    export protfile=`ls ${pathProteinSequences} | grep "${sp}\." | grep fa`

    export speciesList=${sp},${speciesList}
    export pathsFastaProteins=${pathProteinSequences}/${protfile},${pathsFastaProteins}
    export pathsDiamondResults=${pathResults}/${file},${pathsDiamondResults}

    echo ${sp} ${protfile}
done   

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_filter_diamond

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=filter_${target}" >>  ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --output=${pathScripts}/std_output_filter_diamond_${target}.txt" >>  ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --error=${pathScripts}/std_error_filter_diamond_${target}.txt" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --time=1:00:00" >> ${pathScripts}/bsub_script_filter_diamond
fi

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=filter_${refsp}" >>  ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --output=${pathScripts}/std_output_filter_diamond_${refsp}.txt" >>  ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --error=${pathScripts}/std_error_filter_diamond_${refsp}.txt" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --mem=2G" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_filter_diamond
    echo "#SBATCH --time=1:00:00" >> ${pathScripts}/bsub_script_filter_diamond
fi

echo "perl ${pathScripts}/filter.diamond.results.pl --speciesList=${speciesList} --pathsFastaProteins=${pathsFastaProteins} --pathsDiamondResults=${pathsDiamondResults}  --minProteinFraction=0.5 --maxEValue=0.001 --maxGapFraction=0.1 --pathOutput=${pathResults}/SignificantHits_MinProteinFraction0.5_MaxEvalue0.001_MaxGapGraction0.1.txt">> ${pathScripts}/bsub_script_filter_diamond

echo "perl ${pathScripts}/filter.diamond.results.pl --speciesList=${speciesList} --pathsFastaProteins=${pathsFastaProteins} --pathsDiamondResults=${pathsDiamondResults} --minProteinFraction=0.25 --maxEValue=0.001 --maxGapFraction=0.1 --pathOutput=${pathResults}/SignificantHits_MinProteinFraction0.25_MaxEvalue0.001_MaxGapGraction0.1.txt">> ${pathScripts}/bsub_script_filter_diamond

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
     sbatch ${pathScripts}/bsub_script_filter_diamond
fi

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_filter_diamond
    ${pathScripts}/bsub_script_filter_diamond
fi

#########################################################################
