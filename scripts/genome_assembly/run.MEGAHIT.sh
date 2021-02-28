#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathData=${path}/data/WGS
export pathResults=${path}/results/genome_assembly/MEGAHIT
export pathScripts=${path}/scripts/genome_assembly

## megahit: MEGAHIT v1.2.9

#########################################################################

export pathR1=""
export pathR2=""

for file in `ls ${pathData} | grep _R1_001_trimmed.fastq.gz`
do
    export prefix=`basename ${file} _R1_001_trimmed.fastq.gz`
    export pathR1=${pathData}/${prefix}_R1_001_trimmed.fastq.gz,${pathR1}
    export pathR2=${pathData}/${prefix}_R2_001_trimmed.fastq.gz,${pathR2}
done

#########################################################################

echo ${pathR1}
echo ${pathR2}

#########################################################################

if [ ${cluster} = "pbil" ]; then
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --job-name=mh" >>  ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --output=${pathScripts}/std_output_MEGAHIT.txt" >>  ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --error=${pathScripts}/std_error_MEGAHIT.txt" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --cpus-per-task=8" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_megahit

    echo "megahit -1 ${pathR1} -2 ${pathR2} -t 8 --no-mercy --min-count 3 -m 1e10 --out-prefix final -o ${pathResults} --tmp-dir ${pathResults}/tmp" >> ${pathScripts}/bsub_script_megahit
    
    sbatch ${pathScripts}/bsub_script_megahit
fi

#########################################################################

if [ ${cluster} = "cloud" ]; then
    megahit -1 ${pathR1} -2 ${pathR2} -t 8 --no-mercy --min-count 3 -m 1e10 --out-prefix final -o ${pathResults} --tmp-dir ${pathResults}/tmp
fi

#########################################################################
