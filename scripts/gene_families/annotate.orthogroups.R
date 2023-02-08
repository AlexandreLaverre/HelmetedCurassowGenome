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
  
  og=read.table(paste(pathOrthoGroups, node,".tsv",sep=""),h=T, stringsAsFactors=F, sep="\t")

  nbg.ref=unlist(lapply(og[,refsp], function(x) length(unlist(strsplit(x, split=",")))))
  geneid.ref=rep(NA, nrow(og))
  geneid.ref[which(nbg.ref==1)]=unlist(lapply(og[which(nbg.ref==1),refsp], function(x) unlist(strsplit(x, split="\\."))[1]))

  this.names=get(paste("genenames.",refann,sep=""))
  name.ref=this.names[geneid.ref,2]
  
  this.ortho=get(paste("ortho.human.",refann,sep=""))

  geneid.human=rep(NA, nrow(og))
  geneid.human[which(geneid.ref%in%rownames(this.ortho))]=this.ortho[geneid.ref[which(geneid.ref%in%rownames(this.ortho))], "Gene.stable.ID"]
  name.human=rep(NA, nrow(og))
  name.human[which(geneid.human%in%rownames(genenames.human))]=genenames.human[geneid.human[which(geneid.human%in%rownames(genenames.human))],2]

  results=data.frame("HOG"=og$HOG, "OG"=og$OG, "Reference.GeneID"=geneid.ref, "Reference.GeneName"=name.ref, "HumanOrthologue.GeneID"=geneid.human, "HumanOrthologue.GeneName"=name.human)

  write.table(results, paste(pathOrthoGroups, node, "_annotations.tsv",sep=""),row.names=F, col.names=T, sep="\t", quote=F)
  
}

##########################################################################

