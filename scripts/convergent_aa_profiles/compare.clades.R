################################################################

pathResults="../../results/coding_gene_evolution/all_species/"

################################################################

best.birds=read.table(paste(pathResults, "birds/pelican_output_general/best_sites_annotated.tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")

best.squamates=read.table(paste(pathResults, "squamates/pelican_output_general/best_sites_annotated.tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")

################################################################

all.birds=read.table(paste(pathResults, "birds/pelican_output_general/all_sites_annotated.tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")

all.squamates=read.table(paste(pathResults, "squamates/pelican_output_general/all_sites_annotated.tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")

################################################################

signif.birds=best.birds[which(best.birds$aagtr_pval<0.01),]

signif.squamates=best.squamates[which(best.squamates$aagtr_pval<0.01),]

################################################################
