#!/bin/bash

export cluster=$1
export geneset=$2
export dataset=$3

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
    export pathTools=/mnt/mydatalocal/Tools/OrthoFinder/tools
fi


if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
    export pathTools=/beegfs/home/necsulea/Tools/OrthoFinder/tools
fi

export pathCDS=${path}/data/coding_sequences
export pathGeMoMa=${path}/results/genome_annotation
export pathGeneFamilies=${path}/results/gene_families/OrthoFinder/${geneset}/iqtree
export pathResults=${path}/results/coding_gene_evolution/${geneset}/${dataset}
export pathScripts=${path}/scripts/coding_gene_evolution

##########################################################################

if [ -e ${pathResults} ]; then
    echo "output path already there"
else
    mkdir -p ${pathResults}
fi

##########################################################################

export pathsCDS=""
export speciesList=""

for sp in `grep Ensembl ${pathScripts}/species_list.txt | cut -f 1`
do
    export file=`ls ${pathCDS}/Ensembl103 | grep all.fa | grep ${sp}'\.'`

     if [ -e ${pathCDS}/Ensembl103/primary_transcripts/${file} ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathCDS}/Ensembl103/${file}
    fi

    export speciesList=${sp},${speciesList}
    export pathsCDS=${pathCDS}/Ensembl103/primary_transcripts/${file},${pathsCDS}
done

##########################################################################

## add NCBI species

for sp in `grep NCBI ${pathScripts}/species_list.txt | cut -f 1`
do
    if [ -e ${pathGeMoMa}/${sp}/NCBI/GeMoMa/combined/primary_transcripts/combined_annotations_NCBI_GeMoMa_formatted.cds.fa ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathGeMoMa}/${sp}/NCBI/GeMoMa/combined/combined_annotations_NCBI_GeMoMa_formatted.cds.fa
    fi

    export speciesList=${sp},${speciesList}
    export pathsCDS=${pathGeMoMa}/${sp}/NCBI/GeMoMa/combined/primary_transcripts/combined_annotations_NCBI_GeMoMa_formatted.cds.fa,${pathsCDS}
done

##########################################################################

## Basiliscus and Pauxi

for sp in Pauxi_pauxi Basiliscus_vittatus
do
    if [ -e ${pathGeMoMa}/${sp}/MEGAHIT_RAGOUT/GeMoMa/combined/primary_transcripts/filtered_GeMoMa_annotations_formatted.cds.fa ]; then
	echo "primary transcripts already done"
    else
	python ${pathTools}/primary_transcript.py ${pathGeMoMa}/${sp}/MEGAHIT_RAGOUT/GeMoMa/combined/filtered_GeMoMa_annotations_formatted.cds.fa
    fi

    export speciesList=${sp},${speciesList}
    export pathsCDS=${pathGeMoMa}/${sp}/MEGAHIT_RAGOUT/GeMoMa/combined/primary_transcripts/filtered_GeMoMa_annotations_formatted.cds.fa,${pathsCDS}
done

##########################################################################

if [ -e ${pathResults}/CDS ]; then
    echo "dir output exists"
else
    mkdir ${pathResults}/CDS
fi

if [ ${geneset} = "without_chameleons" ]; then
    if [ ${dataset} = "all_species" ]; then
	export pathOrthogroups=`ls ${pathGeneFamilies}/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv`
    fi
    
    if [ ${dataset} = "birds" ]; then
	export pathOrthogroups=`ls ${pathGeneFamilies}/*/Phylogenetic_Hierarchical_Orthogroups/N2.tsv`
    fi
fi


if [ ${geneset} = "all_species" ]; then
    if [ ${dataset} = "all_species" ]; then
	export pathOrthogroups=`ls ${pathGeneFamilies}/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv`
    fi
    
    if [ ${dataset} = "birds" ]; then
	export pathOrthogroups=`ls ${pathGeneFamilies}/*/Phylogenetic_Hierarchical_Orthogroups/N2.tsv`
    fi
fi

echo "path ortho" ${pathOrthogroups}

dos2unix ${pathOrthogroups}

##########################################################################

perl ${pathScripts}/extract.coding.sequences.pl --speciesList=${speciesList} --pathsCDS=${pathsCDS}  --pathOrthogroups=${pathOrthogroups} --requiredSpecies=NA --minNbSpecies=5 --dirOutput=${pathResults}/CDS

##########################################################################
##########################################################################

