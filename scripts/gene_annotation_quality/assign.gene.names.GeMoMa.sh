#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined
export pathAnnot=${path}/data/genome_annotations/Ensembl103
export pathScripts=${path}/scripts/gene_annotation

#########################################################################

export subdirOrtho=`ls ${pathResults}/OrthoFinder_final/OrthoFinder`
export dirOrtho=${pathResults}/OrthoFinder_final/OrthoFinder/${subdirOrtho}/Orthologues/Orthologues_Pauxi_pauxi

## we have to run dos2unix on results

dos2unix ${dirOrtho}/*tsv

#########################################################################

export tgSpecies=Mus_musculus,Gallus_gallus,Anas_platyrhynchos_platyrhynchos

#########################################################################

perl ${pathScripts}/assign.gene.names.pl --refSpecies=Pauxi_pauxi --tgSpecies=${tgSpecies} --dirAnnot=${pathAnnot} --dirOrtho=${dirOrtho} --pathOutput=${pathResults}/homology_inferred_gene_names.txt

#########################################################################
