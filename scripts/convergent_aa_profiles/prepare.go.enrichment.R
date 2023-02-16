################################################################

pathResults="../../results/coding_gene_evolution/all_species/"

################################################################

for(dataset in c("birds", "squamates")){
    print("reading all sites")
    all.sites=read.table(paste(pathResults, dataset, "/pelican_output_general/all_sites_annotated.tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")
    print("done")

    ## select genes that are associated with a single human gene name
    all.sites=all.sites[which(all.sites$HumanGeneName!=""),]
    all.sites=all.sites[grep(",", all.sites$HumanGeneName, invert=T),]

    min.pval.gene=tapply(all.sites$aagtr_pval, as.factor(all.sites$HumanGeneName), min)
    df.min.pval.gene=data.frame("GeneName"=levels(as.factor(all.sites$HumanGeneName)), "MinPValue"=min.pval.gene)
    df.min.pval.gene=df.min.pval.gene[order(df.min.pval.gene$MinPValue),]

    write.table(df.min.pval.gene, file=paste(pathResults, dataset, "/pelican_output_general/minimum_pvalue_per_gene.tsv",sep=""),sep="\t", row.names=F, col.names=T)

    writeLines(df.min.pval.gene$GeneName,con=paste(pathResults, dataset, "/pelican_output_general/genes_sorted_by_minimum_pvalue.txt",sep=""))
}

################################################################
