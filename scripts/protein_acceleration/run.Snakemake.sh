#!/bin/bash

export SpName=$1		# i.e: Duck
export GenomeName=$2		# i.e: Anas_platyrhynchos_platyrhynchos
export Prefix=$3		# i.e: Duck_ATAC_allpeaks
export InterestFile=$4		# i.e: /beegfs/data/alaverre/IPLOSS/results/peaks_evolution/simple_bed/Anas_platyrhynchos_platyrhynchos_ATACSeq_test_1000.txt
export HALFile=$5		# i.e: avian_366/366-avian.hal
export NbPart=$6		# i.e: 10
export minAlign=$7		# i.e: 0.5
export Cluster=$8		# i.e: pbil

##################################################################

if [ ${Cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/IPLOSS
fi


if [ ${Cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
fi

export pathHAL=${path}/results/whole_genome_alignments/${HALFile}
export pathLog=${path}/scripts/conserved_element_evolution/log
mkdir -p ${pathLog}

source ~/.bashrc
conda activate /beegfs/data/alaverre/Tools/envs/envs/snakemake.6

##################################################################

snakemake -j 50 --config SpName=${SpName} RefSp=${GenomeName} Prefix=${Prefix} InterestFile=${InterestFile} pathHAL=${pathHAL} NbPart=${NbPart} Cluster=${Cluster} minAlign=${minAlign} --rerun-incomplete --cluster "sbatch -p normal -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################
