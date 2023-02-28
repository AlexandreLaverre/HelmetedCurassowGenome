###################################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/Species_Tree/"

pathResults="../../results/gene_family_evolution/all_species/birds/"

###################################################################################

library(ape)

###################################################################################

sptree=read.tree(paste(pathOrthoFinder, "SpeciesTree_rooted.txt",sep=""))

###################################################################################

anc.birds=getMRCA(sptree, tip=c("Struthio_camelus_australis", "Gallus_gallus"))

sptree=extract.clade(sptree, node=anc.birds)

write.tree(sptree, file=paste(pathResults, "species_tree.txt",sep=""))

###################################################################################

helmeted.species=c("Numida_meleagris", "Pauxi_pauxi", "Anser_cygnoides", "Anseranas_semipalmata", "Buceros_rhinoceros", "Bucorvus_abyssinicus", "Balearica_regulorum", "Casuarius_casuarius", "Basiliscus_vittatus", "Chamaeleo_calyptratus")

index.tips=which(sptree$tip.label%in%helmeted.species)

edge.length=rep(1, nrow(sptree$edge))
edge.length[which(sptree$edge[,2]%in%index.tips)]=2

index.bucanc=getMRCA(sptree, c("Buceros_rhinoceros", "Bucorvus_abyssinicus"))
edge.length[which(sptree$edge[,2]%in%index.bucanc)]=2

sptree$edge.length=edge.length

###################################################################################

write.tree(sptree, file=paste(pathResults, "species_tree_annotated.txt",sep=""))

###################################################################################

