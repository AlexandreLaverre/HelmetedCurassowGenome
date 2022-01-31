#!/bin/bash

export sample=$1
export cluster=$2
export nthreads=$3

################################################################

if [ ${cluster} = "in2p3" ]; then
    export path="/sps/biometr/necsulea/HelmetedCurassowGenome"
fi

if [ ${cluster} = "pbil" ]; then
    export path="/beegfs/data/necsulea/HelmetedCurassowGenome"
fi
    
export pathScripts=${path}/scripts/FaceBase_data_analysis

################################################################

if [ ${cluster} = "in2p3" ]; then
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_exp
    echo "#SBATCH --job-name=exp_${sample}" >>  ${pathScripts}/bsub_script_exp
    echo "#SBATCH --output=${pathScripts}/std_output_${sample}.txt" >>  ${pathScripts}/bsub_script_exp
    echo "#SBATCH --error=${pathScripts}/std_error_${sample}.txt" >> ${pathScripts}/bsub_script_exp
    echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_exp
    echo "#SBATCH --cpus-per-task=${nthreads}" >> ${pathScripts}/bsub_script_exp
    echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_exp
    
    echo "Rscript --vanilla ${pathScripts}/count.reads.subread.R ${sample} ${nthreads}" >> ${pathScripts}/bsub_script_exp

    sbatch ${pathScripts}/bsub_script_exp
fi

################################################################
