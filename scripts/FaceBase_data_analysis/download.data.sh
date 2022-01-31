#!/bin/bash

export index=$1
export cluster=$2

########################################################################################

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

export pathData=${path}/data/RNASeq/Mouse
export pathScripts=${path}/scripts/FaceBase_data_analysis

########################################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_download

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=d${i}" >>  ${pathScripts}/bsub_script_download
    echo "#SBATCH --output=${pathScripts}/std_output_download_${i}.txt" >>  ${pathScripts}/bsub_script_download
    echo "#SBATCH --error=${pathScripts}/std_error_download_${i}.txt" >> ${pathScripts}/bsub_script_download
    echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_download
    echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_download
    echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_download

    echo "cd ${pathData}" >>  ${pathScripts}/bsub_script_download

    echo -n "wget " >>  ${pathScripts}/bsub_script_download

    sed -n ${i}p ${pathData}/FaceBase_samples.txt >>  ${pathScripts}/bsub_script_download
    
    #sbatch ${pathScripts}/bsub_script_download
fi

########################################################################################
