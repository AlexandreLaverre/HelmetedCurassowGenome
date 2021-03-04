#!/bin/bash

export kmer=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_assembly/SOAPdenovo/
export pathScripts=${path}/scripts/genome_assembly

#########################################################################

# Version 2.04: released on July 13th, 2012
# Compile Mar 22 2020	15:58:14

#########################################################################

if [ ${cluster} = "pbil" ]; then
    
    echo "#!/bin/bash " > ${pathScripts}/bsub_script_test_SOAPdenovo
    
    echo "#SBATCH --job-name=SOAP_${kmer}" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    echo "#SBATCH --output=${pathScripts}/std_out_SOAP_${kmer}" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    echo "#SBATCH --error=${pathScripts}/std_err_SOAP_${kmer}" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    echo "#SBATCH --cpus-per-task=8" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    echo "#SBATCH --mem=20G" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    
    echo "soapdenovo2-127mer all -s ${pathScripts}/configFile_SOAPdenovo_${cluster} -o ${pathResults}/kmer${kmer} -p 8 -a 20 -K ${kmer}" >>  ${pathScripts}/bsub_script_test_SOAPdenovo
    
    sbatch ${pathScripts}/bsub_script_test_SOAPdenovo
    
fi

#########################################################################

if [ ${cluster} = "cloud" ]; then
    soapdenovo2-127mer all -s ${pathScripts}/configFile_SOAPdenovo_${cluster} -o ${pathResults}/kmer${kmer} -p 12 -a 120 -K ${kmer}
fi

#########################################################################
