#!/bin/bash

export cluster=$1

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
export pathAnnot=${path}/results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined
export pathResults=${path}/results/coding_gene_evolution
export pathScripts=${path}/scripts/coding_gene_evolution

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

export speciesList=Penelope_pileata,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013396635.1_ASM1339663v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Alectura_lathami,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013399715.1_ASM1339971v1_cds_from_genomic.fna.gz,${pathsCDS} 

##########################################################################

export speciesList=Casuarius_casuarius,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013396415.1_ASM1339641v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Anseranas_semipalmata,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013399115.1_ASM1339911v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Balearica_regulorum,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCF_000709895.1_ASM70989v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Grus_americana,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013390085.1_ASM1339008v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Bucorvus_abyssinicus,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013398885.1_ASM1339888v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Buceros_rhinoceros,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_000710305.1_ASM71030v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Upupa_epops,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013397515.1_ASM1339751v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

export speciesList=Rhinopomastus_cyanomelas,${speciesList}
export pathsCDS=${pathCDS}/NCBI/GCA_013400115.1_ASM1340011v1_cds_from_genomic.fna.gz,${pathsCDS}

##########################################################################

## Pauxi pauxi

if [ -e ${pathAnnot}/primary_transcripts/final_annotations_formatted.cds.fa ]; then
    echo "primary transcripts already there for Pauxi pauxi"
else
    python ${pathTools}/primary_transcript.py ${pathAnnot}/final_annotations_formatted.cds.fa
fi

export speciesList=Pauxi_pauxi,${speciesList}
export pathsCDS=${pathAnnot}/primary_transcripts/final_annotations_formatted.cds.fa,${pathsCDS}

##########################################################################
##########################################################################

if [ -e ${pathResults}/CDS ]; then
    echo "dir output exists"
else
    mkdir ${pathResults}/CDS
fi

export pathOrthogroups=`ls ${pathResults}/OrthoFinder_iqtree/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv`

echo "path ortho" ${pathOrthogroups}

dos2unix ${pathOrthogroups}

##########################################################################

perl ${pathScripts}/extract.coding.sequences.pl --speciesList=${speciesList} --pathsCDS=${pathsCDS}  --pathOrthogroups=${pathOrthogroups} --requiredSpecies=Pauxi_pauxi --minNbSpecies=5 --dirOutput=${pathResults}/CDS

##########################################################################
##########################################################################

