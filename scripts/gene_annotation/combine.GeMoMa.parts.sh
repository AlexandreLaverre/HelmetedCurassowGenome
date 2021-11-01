#!/bin/bash

export assembly=$1
export ref=$2
export source=$3
export cluster=$4
export threads=$5
export constraint=$6

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
    export pathTools=/beegfs/home/${USER}/Tools
    export version=1.8
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
    export pathTools=/sps/biometr/necsulea/Tools
    export version=1.8
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools
    export version=1.8
fi

export pathGenomeAssembly=${path}/results/genome_assembly/${assembly}
export pathGenomes=${path}/data/genome_sequences/${source}
export pathAnnotations=${path}/data/genome_annotations/${source}
export pathResults=${path}/results/genome_annotation/${assembly}/GeMoMa/${ref}
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

if [ -e ${pathResults} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}
fi


#########################################################################

if [ -e ${pathResults}/final_annotation.gff ]; then
    echo "combined annotations exist, not doing anything"
    exit
fi

#########################################################################

for part in `ls ${pathResults} | grep -v final`
do
    export part=`echo ${annotfile} | cut -f 2 -d '.' `
    
    if [ -e ${pathResults}/${part}/final_annotation.gff ]; then
	echo ${part}
	
	if [ -e ${pathResults}/final_annotation.gff ]; then
	    cat ${pathResults}/${part}/final_annotation.gff | grep -v SOFTWARE | grep -v gff-version >> ${pathResults}/final_annotation.gff
	else
	    cp ${pathResults}/${part}/final_annotation.gff ${pathResults}/final_annotation.gff
	fi
    else
	export nbcds=`grep -c CDS ${pathAnnotations}/parts/${ref}.${part}.gff`

	if [ ${nbcds} = 0 ]; then
	    echo "cannot find results for "${part}", but there are no CDS"
	else
	    echo "cannot find results for "${part}", there are "${nbcds}" CDS"
	    exit
	fi
    fi
done

#########################################################################
