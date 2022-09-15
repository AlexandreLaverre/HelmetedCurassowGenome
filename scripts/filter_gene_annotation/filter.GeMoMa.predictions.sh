#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/HelmetedCurassowGenome
fi

export pathResults=${path}/results/genome_annotation/${sp}/${assembly}/GeMoMa/combined
export pathScripts=${path}/scripts/filter_gene_annotation

#########################################################################

if [ ${assembly} = "MEGAHIT" ]; then
    export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
    export pathAssembly=${pathGenomeAssembly}/final.contigs.fa
fi

#########################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export pathGenomeAssembly=${path}/results/genome_assembly/${sp}/${assembly}
    export pathAssembly=${pathGenomeAssembly}/genome_sequence_renamed_sm.fa
fi

#########################################################################

if [ ${assembly} = "NCBI" ]; then
    export pathGenomeSequences=${path}/data/genome_sequences/${assembly}
    export pathAssembly=${pathGenomeSequences}/${sp}.fa
fi

#########################################################################

export pathDiamondResults=${pathResults}/diamond_results/SignificantHits_MinProteinFraction0.25_MaxEvalue0.001_MaxGapGraction0.1.txt

#########################################################################

## quantile(chicken protein lengths, p=0.01) = 68

for minlen in 70 100
do
    for maxrep in 0.25 0.5
    do

	perl ${pathScripts}/filter.GeMoMa.predictions.pl --pathAnnotGTF=${pathResults}/filtered_predictions.gtf --pathProteins=${pathResults}/filtered_predictions.faa --pathDiamondResults=${pathDiamondResults} --minProteinLength=${minlen} --pathOverlapRepeats=${pathResults}/overlap_repeats.txt --maxFractionRepeats=${maxrep} --source=GeMoMa --pathOutputGTF=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.gtf --pathOutputFasta=${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.faa

	gffread -S -x ${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.cds.fa -g ${pathAssembly} ${pathResults}/filtered_predictions_minDiamondProteinFraction0.25_minLength${minlen}_maxFractionRepeats${maxrep}.gtf
	
    done
done
#########################################################################
