########################################################################

library(ape)

########################################################################

pathResults="../../results/coding_gene_evolution/"

########################################################################

tree=read.tree(paste(pathResults, "OrthoFinder_fasttree/Results_Nov16/Species_Tree/SpeciesTree_rooted.txt",sep=""))

rooted=root(tree, outgroup=c("Struthio_camelus_australis", "Casuarius_casuarius", "Dromaius_novaehollandiae"))

write.tree(rooted, file=paste(pathResults, "species_tree.txt",sep=""))

########################################################################
