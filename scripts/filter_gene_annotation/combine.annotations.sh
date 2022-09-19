#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3

#############################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathAnnot=${path}/data/genome_annotations/${assembly}
export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/filter_gene_annotation

#############################################################################

if [ -e ${pathAnnot}/${sp}.gtf ]; then
    echo "gtf file already there"
else
    gffread -E ${pathAnnot}/${sp}.gff -T -o ${pathAnnot}/${sp}.gtf
fi

#############################################################################

perl ${pathScripts}/select.new.transcripts.pl --pathAnnot1=${pathAnnot}/${sp}.gtf --pathAnnot2=${pathResults}/filtered_GeMoMa_annotations.gtf --chrList=NA --pathOutputDecision=${pathResults}/selected_transcripts_vs_${assembly}.txt

#############################################################################

perl ${pathScripts}/combine.annotations.pl --pathAnnot1=${pathAnnot}/${sp}.gtf --pathAnnot2=${pathResults}/filtered_GeMoMa_annotations.gtf  --chrList=NA --pathSelectedTranscripts=${pathResults}/selected_transcripts_vs_${assembly}.txt --pathOutputGTF=${pathResults}/combined_annotations_${assembly}_GeMoMa.gtf

#############################################################################
