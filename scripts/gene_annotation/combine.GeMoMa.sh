#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
    export pathTools=/beegfs/home/${USER}/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools
fi

export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/GeMoMa

#########################################################################

if [ -e ${pathResults}/combined ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}/combined
fi

#########################################################################

echo "#!/bin/bash " > script_combine_GeMoMa
echo -n "java -jar ${pathTools}/GeMoMa/GeMoMa-1.8.jar CLI GAF " >> script_combine_GeMoMa

#########################################################################

for sp in `ls ${pathResults} | grep -v combined`
do
    if [ -e ${pathResults}/${sp}/final_annotation.gff ]; then
	echo -n "p=${sp} g=${pathResults}/${sp}/final_annotation.gff ">> script_combine_GeMoMa
    fi
done

#########################################################################

echo  "outdir=${pathResults}/combined " >> script_combine_GeMoMa

#########################################################################

chmod a+x script_combine_GeMoMa

./script_combine_GeMoMa

#########################################################################
