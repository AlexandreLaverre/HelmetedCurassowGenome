################################################################
## adapted from original scripts by Louis Duchemin

if(!require("gtfisher")){
    devtools::install_github("https://github.com/lsdch/gtfisher")
}

if(!require("progressr")){
    install.packages("progressr")
}

library(gtfisher)
library(progressr)

################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/Phylogenetic_Hierarchical_Orthogroups/"

og=list("birds"="N2", "squamates"="N1")

for(dataset in c("birds", "squamates")){

    this.og=og[[dataset]]
    og.annot=read.table(paste(pathOrthoFinder, this.og, "_annotations.tsv",sep=""),h=T, stringsAsFactors=F,sep="\t", quote="\"")
    rownames(og.annot)=paste(og.annot$HOG, og.annot$OG, sep="_")

    pathResults=paste("../../results/coding_gene_evolution/all_species/", dataset, "/pelican_output_general/", sep="")

    ## Load site-level predictions for all alignments
    site_pvals = readr::read_tsv(paste(pathResults, "all_sites.tsv",sep=""))

    ## Provide loaded dataframe, and columns for alignment and site p-values, while filtering out constant sites (naa == 1).

    ## with progress bar
    handlers("progress")
    with_progress({
        gene_pvals = gtf_predict(site_pvals, alignment, aagtr_pval, k=3, n=1000, naa > 1)
    })

    gene_pvals=as.data.frame(gene_pvals)

    gene_pvals=cbind(gene_pvals, og.annot[as.character(gene_pvals$alignment),c("Reference.GeneID", "Reference.GeneName", "HumanOrthologue.GeneID", "HumanOrthologue.GeneName")])

    best_sites=read.table(paste(pathResults, "all_sites.tsv",sep=""),h=T, stringsAsFactors=F)
    best_sites=cbind(best_sites, og.annot[as.character(best_sites$alignment),c("Reference.GeneID", "Reference.GeneName", "HumanOrthologue.GeneID", "HumanOrthologue.GeneName")])
}

################################################################
