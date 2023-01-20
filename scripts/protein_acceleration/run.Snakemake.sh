#!/bin/bash

export Prefix=$1		# i.e: all_species
export NbPart=$2		# i.e: 10
export minAlign=$3		# i.e: 0.9
export Cluster=$4		# i.e: pbil

##################################################################

if [ ${Cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

if [ ${Cluster} = "local" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

export pathLog=${path}/scripts/protein_acceleration/log
mkdir -p ${pathLog}

source ~/.bashrc
conda activate /beegfs/data/alaverre/Tools/envs/envs/snakemake.6

##################################################################

snakemake -j 100 --config Prefix=${Prefix} NbPart=${NbPart} Cluster=${Cluster} minAlign=${minAlign} --rerun-incomplete --cluster "sbatch -p normal -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################
