###################################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/Phylogenetic_Hierarchical_Orthogroups/"

pathResults="../../results/gene_family_evolution/all_species/birds/"

###################################################################################

og=read.table(paste(pathOrthoFinder, "N2.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(og)=paste(og$HOG, og$OG, sep="_")

og.annot=read.table(paste(pathOrthoFinder, "/N2_annotations.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(og.annot)=paste(og.annot$HOG, og.annot$OG, sep="_")

###################################################################################

splist=colnames(og)[4:ncol(og)]

nbgenes=t(apply(og[,splist], 1, function(x) unlist(lapply(x, function(y) length(unlist(strsplit(y, split=",")))))))

###################################################################################

cafe5.input=data.frame("Desc"=og.annot[rownames(og), "HumanOrthologue.GeneName"], "Family ID"=rownames(og), nbgenes, check.names=F)
cafe5.input$Desc[which(cafe5.input$Desc=="")]="NA"

###################################################################################

write.table(cafe5.input, file=paste(pathResults, "gene_families.txt",sep=""), row.names=F,col.names=T, sep="\t", quote=F)

###################################################################################
