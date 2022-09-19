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
export pathAnnotScripts=${path}/scripts/gene_annotation
export pathScripts=${path}/scripts/filter_gene_annotation

############################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
fi

#############################################################################

if [ ${assembly} = "NCBI" ]; then
    export pathGenomeSequences=${path}/data/genome_sequences/${assembly}
    export pathAssembly=${pathGenomeSequences}/${sp}.fa
fi

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
## extract fasta files

gffread -S -x ${pathResults}/combined_annotations_${assembly}_GeMoMa.cds.fa -g ${pathAssembly} ${pathResults}/combined_annotations_${assembly}_GeMoMa.gtf

gffread -S -y ${pathResults}/combined_annotations_${assembly}_GeMoMa.faa -g ${pathAssembly} ${pathResults}/combined_annotations_${assembly}_GeMoMa.gtf

#############################################################################

## add gene id and transcript id in fasta file

perl ${pathAnnotScripts}/format.GeMoMa.sequences.pl --pathAnnotGTF=${pathResults}/combined_annotations_${assembly}_GeMoMa.gtf --pathSequences=${pathResults}/combined_annotations_${assembly}_GeMoMa.cds.fa --pathOutput=${pathResults}/combined_annotations_${assembly}_GeMoMa_formatted.cds.fa

perl ${pathAnnotScripts}/format.GeMoMa.sequences.pl --pathAnnotGTF=${pathResults}/combined_annotations_${assembly}_GeMoMa.gtf --pathSequences=${pathResults}/combined_annotations_${assembly}_GeMoMa.faa --pathOutput=${pathResults}/combined_annotations_${assembly}_GeMoMa_formatted.faa

#############################################################################
