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
full.tree$node.label <- NULL

write.tree(full.tree, paste(pathResults, "/species_tree_unrooted_nobootstrap.txt",sep=""))

###########################################################################

helmeted=c("Anseranas_semipalmata", "Numida_meleagris", "Casuarius_casuarius", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

all.species=full.tree$tip.label

phenotype.data=data.frame("species"=all.species, "trait"=rep(NA, length(all.species)))
phenotype.data$trait[which(phenotype.data$species%in%helmeted)]="helmeted"
phenotype.data$trait[which(!(phenotype.data$species%in%helmeted))]="non-helmeted"

write.table(phenotype.data, file=paste(pathResults, "/phenotype_data.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

###########################################################################

names=read.table(paste(pathAnnot, "homology_inferred_gene_names.txt",sep=""), h=T, stringsAsFactors=F)
rownames(names)=names$GeneID

###########################################################################

N0=read.table(paste(pathOrthoGroups, "N0.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")

###########################################################################

for(i in 1:nrow(N0)){
  hog=N0$HOG[i]
  og=N0$OG[i]

  path=paste(pathCDS, hog,"_",og,".aln.best.fas",sep="")

  if(file.exists(path)){
    id.pauxi=N0[i,"Pauxi_pauxi"]

    if(id.pauxi%in%rownames(names)){
      name.pauxi=names[id.pauxi, "GeneName.Pauxi_pauxi"]

      if(!is.na(name.pauxi)){
        system(paste("cp ",path, " ",pathOutput, "/",id.pauxi,"_",name.pauxi,".fa",sep=""))
      } else{
        system(paste("cp ",path, " ",pathOutput, "/",id.pauxi,".fa",sep=""))
      }
    } else{
      system(paste("cp ",path, " ",pathOutput, "/",id.pauxi,".fa",sep=""))
    }
  }

  if(i%%1000==0){
    print(i)
  }
}

###########################################################################



