#########################################################################################################

pathAnnot="../../data/ensembl_annotations/Human/"
pathPelican="../../results/coding_gene_evolution/all_species/"
pathLiterature="../../data/data_from_literature/"


dataset="birds"

minseq=5
minaa=2
maxpval=0.01

#########################################################################################################

cc=system(paste("grep cellular_component$ ", pathAnnot, "GeneOntology_Ensembl103.txt",sep=""), intern=T)

cc.gene=unlist(lapply(cc, function(x) unlist(strsplit(x, split="\t"))[1]))
cc.cat=unlist(lapply(cc, function(x) unlist(strsplit(x, split="\t"))[3]))

#########################################################################################################

all.sites=read.table(paste(pathPelican, dataset, "/pelican_output_general/all_sites_annotated.tsv",sep=""), h=T, stringsAsFactors=F,sep="\t", quote="\"") ##  20582600 sites

all.sites=all.sites[which(all.sites$nseq>=minseq & all.sites$naa>=minaa),] ##  3155736 sites

all.sites=all.sites[which(all.sites$HumanGeneID!=""),]
all.sites=all.sites[grep(",", all.sites$HumanGeneID, invert=T),]
all.sites=all.sites[which(all.sites$HumanGeneID%in%cc.gene),]

nb.sites.per.gene=as.numeric(table(all.sites$HumanGeneID))
names(nb.sites.per.gene)=levels(as.factor(all.sites$HumanGeneID))

#########################################################################################################

all.sites=all.sites[order(all.sites$aagtr_pval),]

#########################################################################################################

signif.gene=unique(all.sites$HumanGeneID[which(all.sites$aagtr_pval<maxpval)])
notsignif.gene=unique(all.sites$HumanGeneID[which(all.sites$aagtr_pval>=maxpval)])

nb.signif.sites=length(which(all.sites$aagtr_pval<maxpval))

#########################################################################################################

enrichment=list()

for(cat in c("extracellular matrix", "centriole")){

    this.cc.gene=cc.gene[which(cc.cat==cat)]

    prop.signif.cat=length(intersect(this.cc.gene, signif.gene))/length(signif.gene)
    prop.notsignif.cat=length(intersect(this.cc.gene, notsignif.gene))/length(notsignif.gene)

    enrichment[[cat]]=list("prop.signif"=prop.signif.cat,"prop.notsignif"=prop.notsignif.cat)

    rand.prop.signif=c()

    for(i in 1:1000){
        print(i)

        fake.signif.sites=sample(1:nrow(all.sites),size=nb.signif.sites, re=FALSE)
        fake.signif.gene=unique(all.sites[fake.signif.sites, "HumanGeneID"])

        fake.prop.signif.cat=length(intersect(this.cc.gene, fake.signif.gene))/length(fake.signif.gene)

        rand.prop.signif=c(rand.prop.signif, fake.prop.signif.cat)
    }

    enrichment[[cat]][["rand.prop.signif"]]=rand.prop.signif
}

#########################################################################################################

signif.name=unique(all.sites$HumanGeneName[which(all.sites$aagtr_pval<maxpval)])
notsignif.name=unique(all.sites$HumanGeneName[which(all.sites$aagtr_pval>=maxpval)])
signif.name=setdiff(signif.name, "")
notsignif.name=setdiff(notsignif.name, "")

cil.genes=readLines(paste(pathLiterature, "ciliopathy_genes_Reiter_Leroux_2017.txt", sep=""))

prop.signif.cil=length(intersect(cil.genes, signif.name))/length(signif.name)
prop.notsignif.cil=length(intersect(cil.genes, notsignif.name))/length(notsignif.name)

rand.prop.signif=c()

for(i in 1:1000){
    print(i)

    fake.signif.sites=sample(1:nrow(all.sites),size=nb.signif.sites, re=FALSE)
    fake.signif.name=unique(all.sites[fake.signif.sites, "HumanGeneName"])

    fake.prop.signif.cil=length(intersect(cil.genes, fake.signif.name))/length(fake.signif.name)

    rand.prop.signif=c(rand.prop.signif, fake.prop.signif.cil)
}

enrichment[["ciliopathy"]]=list("prop.signif"=prop.signif.cil, "prop.notsignif"=prop.notsignif.cil, "rand.prop.signif"=rand.prop.signif)

#########################################################################################################

save(enrichment, file=paste(pathPelican, dataset, "/pelican_output_general/functional.enrichment.RData",sep=""))

#########################################################################################################

