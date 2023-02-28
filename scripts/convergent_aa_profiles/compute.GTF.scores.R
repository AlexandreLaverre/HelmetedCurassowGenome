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

for(dataset in c("birds", "squamates")){

    pathResults=paste("../../results/coding_gene_evolution/all_species/", dataset, "/pelican_output_general/", sep="")

    ## Load site-level predictions for all alignments
    site.pvals = readr::read_tsv(paste(pathResults, "all_sites_annotated.tsv",sep=""))

    site.info=as.data.frame(site.pvals)

    aln.info=data.frame("alignment"=site.info$alignment, "name"=site.info$HumanGeneName)
    aln.info=aln.info[which(!duplicated(aln.info$alignment)),]
    rownames(aln.info)=aln.info$alignment

    ## Provide loaded dataframe, and columns for alignment and site p-values, while filtering out constant sites (naa == 1).

    ## with progress bar
    handlers("progress")
    with_progress({
        gene.pvals = gtf.predict(site.pvals, alignment, aagtr.pval, k=1, n=1000, naa > 1, nseq > 4)
    })

    gene.pvals=as.data.frame(gene.pvals)
    gene.pvals$HumanGeneName=aln.info[gene.pvals$alignment, "name"]

}

################################################################
