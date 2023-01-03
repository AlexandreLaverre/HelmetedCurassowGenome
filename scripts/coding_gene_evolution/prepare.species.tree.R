########################################################################

library(ape)

########################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/iqtree/Results_Sep21/"
pathResults="../../results/coding_gene_evolution/"

########################################################################

tree=read.tree(paste(pathOrthoFinder, "/Species_Tree/SpeciesTree_rooted.txt",sep=""))

rooted=root(tree, outgroup=c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis"))

write.tree(rooted, file=paste(pathResults, "species_tree.txt",sep=""))

########################################################################
