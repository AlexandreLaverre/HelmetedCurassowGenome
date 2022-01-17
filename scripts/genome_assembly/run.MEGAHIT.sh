#!/bin/bash

export species=$1
export cluster=$2
export nthreads=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathData=${path}/data/WGS/${species}
export pathResults=${path}/results/genome_assembly/${species}
export pathScripts=${path}/scripts/genome_assembly

## megahit: MEGAHIT v1.2.9

#########################################################################

export pathR1=""
export pathR2=""

for file in `ls ${pathData} | grep R1 | grep trimmed.fastq.gz`
do
    export pathR1=${pathData}/${file},${pathR1}
done

for file in `ls ${pathData} | grep R2 | grep trimmed.fastq.gz`
do
    export pathR2=${pathData}/${file},${pathR2}
done

#########################################################################

echo ${pathR1}
echo ${pathR2}

#########################################################################

if [ -e ${pathResults} ]; then
    echo "output dir already there"
else
    mkdir -p ${pathResults}
fi

if [ -e ${pathResults}/tmp ]; then
    echo "tmp dir already there"
else
    mkdir -p ${pathResults}/tmp
fi

########################################################################

if [ ${cluster} = "pbil" ]; then
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --job-name=mh" >>  ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --output=${pathScripts}/std_output_MEGAHIT.txt" >>  ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --error=${pathScripts}/std_error_MEGAHIT.txt" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --cpus-per-task=${nthreads}" >> ${pathScripts}/bsub_script_megahit
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_megahit

    echo "megahit -1 ${pathR1} -2 ${pathR2} -t ${nthreads} --no-mercy --min-count 3 -m 1e10 --out-prefix final -o ${pathResults}/MEGAHIT --tmp-dir ${pathResults}/tmp" >> ${pathScripts}/bsub_script_megahit

    sbatch ${pathScripts}/bsub_script_megahit
fi

#########################################################################

if [ ${cluster} = "cloud" ]; then
    megahit -1 ${pathR1} -2 ${pathR2} -t ${nthreads} --no-mercy --min-count 3 -m 1e10 --out-prefix final -o ${pathResults}/MEGAHIT --tmp-dir ${pathResults}/tmp
fi

#########################################################################
