###########################################################################

pathResults="../../results/coding_gene_evolution"
pathOrthoGroups=paste(pathResults, "/OrthoFinder_iqtree/Results_Nov17/Phylogenetic_Hierarchical_Orthogroups/",sep="")
pathCDS=paste(pathResults, "/CDS/",sep="")
pathOutput=paste(pathResults, "/data_for_pelican/",sep="")
pathAnnot="../../results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined/"

###########################################################################

library(ape)

###########################################################################

full.tree=read.tree(paste(pathResults, "/species_tree_rooted.txt",sep=""))

if(file.exists(paste(pathResults, "/species_tree_nobootstrap.txt",sep=""))){
  print("tree already there")
}  else {
  full.tree$node.label <- NULL
  write.tree(full.tree, paste(pathResults, "/species_tree_nobootstrap.txt",sep=""))
}

###########################################################################

withIP=c("Struthio_camelus_australis", "Dromaius_novaehollandiae", "Casuarius_casuarius", "Anser_brachyrhynchus", "Anseranas_semipalmata","Anas_platyrhynchos_platyrhynchos", "Anser_cygnoides", "Penelope_pileata", "Pauxi_pauxi")


all.species=full.tree$tip.label

phenotype.data=data.frame("species"=all.species, "trait"=rep(NA, length(all.species)))
phenotype.data$trait[which(!(phenotype.data$species%in%withIP))]="without_IP"

phenotype.data$trait[which(phenotype.data$species%in%withIP)]="with_IP"

write.table(phenotype.data, file=paste(pathResults, "/phenotype_data_IP.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)


###########################################################################
