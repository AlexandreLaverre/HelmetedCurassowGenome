###########################################################################

pathResults="../../results/coding_gene_evolution"
pathOrthoGroups=paste(pathResults, "/OrthoFinder_iqtree/Results_Nov17/Phylogenetic_Hierarchical_Orthogroups/",sep="")
pathCDS=paste(pathResults, "/CDS/",sep="")
pathOutput=paste(pathResults, "/data_for_pelican/",sep="")
pathAnnot="../../results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined/"

###########################################################################

library(ape)

###########################################################################

if(file.exists(paste(pathResults, "/species_tree_nobootstrap.txt",sep=""))){
  print("tree already there")
}  else {
  full.tree=read.tree(paste(pathResults, "/species_tree_rooted.txt",sep=""))
  full.tree$node.label <- NULL
  
  write.tree(full.tree, paste(pathResults, "/species_tree_nobootstrap.txt",sep=""))
}

###########################################################################

for(dataset in c("all_species", "by_category")){

  helmeted=c("Anseranas_semipalmata", "Numida_meleagris", "Casuarius_casuarius", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

  upper_beak=c("Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")
  dorsal_neurocranium=c("Numida_meleagris", "Casuarius_casuarius")
  frontal_area=c("Anseranas_semipalmata", "Balearica_regulorum")

  all.species=full.tree$tip.label
  
  phenotype.data=data.frame("species"=all.species, "trait"=rep(NA, length(all.species)))
  phenotype.data$trait[which(!(phenotype.data$species%in%helmeted))]="non-helmeted"

  if(dataset == "all_species"){
    phenotype.data$trait[which(phenotype.data$species%in%helmeted)]="helmeted"
  }

  if(dataset == "by_category"){
    phenotype.data$trait[which(phenotype.data$species%in%upper_beak)]="upper_beak"
    phenotype.data$trait[which(phenotype.data$species%in%dorsal_neurocranium)]="dorsal_neurocranium"
    phenotype.data$trait[which(phenotype.data$species%in%frontal_area)]="frontal_area"
  }
  
  write.table(phenotype.data, file=paste(pathResults, "/phenotype_data_",dataset,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

###########################################################################
