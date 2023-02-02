###################################################################

library(ape)

###################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/"
pathFigures="../../results/figures/"

###################################################################

sptree=read.tree(file=paste(pathOrthoFinder, "Species_Tree/SpeciesTree_rooted.txt",sep=""))

###################################################################

save(sptree, file=paste(pathFigures, "RData/species.phylogeny.RData",sep=""))

###################################################################
