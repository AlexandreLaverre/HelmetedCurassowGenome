##########################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/"
pathOrthoGroups=paste(pathOrthoFinder, "Phylogenetic_Hierarchical_Orthogroups/", sep="")
pathAnnot="../../data/ensembl_annotations/"
pathEnsemblHomology="../../data/ensembl_homology/"

##########################################################################

nodes=list("birds"="N2", "squamates"="N1")
refspecies=list("birds"="Gallus_gallus", "squamates"="Anolis_carolinensis")
refannot=list("birds"="chicken", "squamates"="lizard")

##########################################################################

genenames.chicken=read.table(paste(pathAnnot, "Chicken/GeneNames_Ensembl103.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(genenames.chicken)=genenames.chicken[,1]

genenames.lizard=read.table(paste(pathAnnot, "Anolis_carolinensis/GeneNames_Ensembl103.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(genenames.lizard)=genenames.lizard[,1]

genenames.human=read.table(paste(pathAnnot, "Human/GeneNames_Ensembl103.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(genenames.human)=genenames.human[,1]

##########################################################################

ortho.human.chicken=read.table(paste(pathEnsemblHomology, "Homologues_Human_Chicken_Ensembl103.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
ortho.human.chicken=ortho.human.chicken[which(ortho.human.chicken[,3]=="ortholog_one2one"),]
rownames(ortho.human.chicken)=ortho.human.chicken[, "Chicken.gene.stable.ID"]

ortho.human.lizard=read.table(paste(pathEnsemblHomology, "Homologues_Human_Lizard_Ensembl103.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
ortho.human.lizard=ortho.human.lizard[which(ortho.human.lizard[,3]=="ortholog_one2one"),]
rownames(ortho.human.lizard)=ortho.human.lizard[, "Anole.lizard.gene.stable.ID"]

##########################################################################

for(splist in c("birds", "squamates")){
  node=nodes[[splist]]
  refsp=refspecies[[splist]]
  refann=refannot[[splist]]

  this.names=get(paste("genenames.",refann,sep=""))
  this.ortho=get(paste("ortho.human.",refann,sep=""))

  og=read.table(paste(pathOrthoGroups, node,".tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")

  list.geneid.ref=lapply(og[,refsp], function(x) unlist(strsplit(x, split=", ")))
  list.geneid.ref=lapply(list.geneid.ref, function(x) unlist(lapply(x, function(y) unlist(strsplit(y, split="\\."))[1])))
  list.genename.ref=lapply(list.geneid.ref, function(x) this.names[intersect(x, rownames(this.names)),2])

  geneid.ref=unlist(lapply(list.geneid.ref, function(x) paste(x, collapse=",")))
  genename.ref=unlist(lapply(list.genename.ref, function(x) paste(setdiff(x,""), collapse=",")))

  list.geneid.human=lapply(list.geneid.ref, function(x) this.ortho[intersect(x, rownames(this.ortho)), "Gene.stable.ID"])
  list.genename.human=lapply(list.geneid.human,function(x) genenames.human[intersect(x, rownames(genenames.human)),2])

  geneid.human=unlist(lapply(list.geneid.human, function(x) paste(x, collapse=",")))
  genename.human=unlist(lapply(list.genename.human, function(x) paste(setdiff(x,""), collapse=",")))

  results=data.frame("HOG"=og$HOG, "OG"=og$OG, "Reference.GeneID"=geneid.ref, "Reference.GeneName"=genename.ref, "HumanOrthologue.GeneID"=geneid.human, "HumanOrthologue.GeneName"=genename.human)

  write.table(results, paste(pathOrthoGroups, node, "_annotations.tsv",sep=""),row.names=F, col.names=T, sep="\t", quote=F)

}

##########################################################################

