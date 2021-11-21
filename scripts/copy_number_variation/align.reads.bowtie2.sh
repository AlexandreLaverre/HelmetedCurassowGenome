#!/bin/bash

export cluster=$1
export threads=$2

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

##############################################################

export pathData=${path}/data/WGS
export pathResults=${path}/results/bowtie2_alignments
export pathScripts=${path}/scripts/copy_number_variation

export pathBowtieIndex=${path}/data/genome_indexes/MEGAHIT_RAGOUT/genome_sequence_renamed

## bowtie2 2.4.4

##############################################################

## check files

export filesR1=""
export filesR2=""

for file in `ls ${pathData} | grep R1`
do
    export prefix=`basename ${file} _R1_001_trimmed.fastq.gz`
    export filesR1=${pathData}/${prefix}_R1_001_trimmed.fastq.gz,${filesR1}
    export filesR2=${pathData}/${prefix}_R2_001_trimmed.fastq.gz,${filesR2}
done

echo ${filesR1}
echo ${filesR2}

##############################################################

if [ -e ${pathResults}/accepted_hits_allsamples.sam ]||[ -e ${pathResults}/accepted_hits_allsamples.bam ]; then
    echo "already done"
else

    echo "#!/bin/bash " > ${pathScripts}/log/bsub_script_bowtie2

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=bowtie2" >>  ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --partition=normal" >>  ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --output=${pathScripts}/log/std_out_bowtie2" >>  ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --error=${pathScripts}/log/std_err_bowtie2" >>  ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --cpus-per-task=${threads}" >>  ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --time=2:00:00" >>  ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --mem=5G" >>  ${pathScripts}/log/bsub_script_bowtie2
    fi

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=bowtie2" >> ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --output=${pathScripts}/log/std_out_bowtie2" >> ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --error=${pathScripts}/log/std_err_bowtie2" >> ${pathScripts}/log/bsub_script_bowtie2
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/log/bsub_script_bowtie2
    fi
    
    echo "bowtie2 -p ${threads} --no-unal --very-sensitive-local -x ${pathBowtieIndex} -1 ${filesR1} -2 ${filesR2} -S ${pathResults}/accepted_hits_allsamples.sam > ${pathResults}/bowtie2_output 2>&1" >> ${pathScripts}/log/bsub_script_bowtie2

    echo "samtools sort -m 12G -@ ${threads} -o ${pathResults}/accepted_hits_allsamples.bam ${pathResults}/accepted_hits_allsamples.sam " >> ${pathScripts}/log/bsub_script_bowtie2

    if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
	sbatch ${pathScripts}/log/bsub_script_bowtie2
    fi

    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/log/bsub_script_bowtie2
	${pathScripts}/log/bsub_script_bowtie2
    fi
fi


##############################################################
