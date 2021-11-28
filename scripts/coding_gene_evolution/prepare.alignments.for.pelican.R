###########################################################################

pathResults="../../results/coding_gene_evolution"
pathOrthoGroups=paste(pathResults, "/OrthoFinder_iqtree/Results_Nov17/Phylogenetic_Hierarchical_Orthogroups/",sep="")
pathCDS=paste(pathResults, "/CDS/",sep="")
pathOutput=paste(pathResults, "/data_for_pelican/",sep="")
pathAnnot="../../results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined/"

###########################################################################

library(stringr)

###########################################################################

names=read.table(paste(pathAnnot, "homology_inferred_gene_names.txt",sep=""), h=T, stringsAsFactors=F)
rownames(names)=names$GeneID

###########################################################################

N0=read.table(paste(pathOrthoGroups, "N0.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")

###########################################################################

for(i in 1:nrow(N0)){
  hog=N0$HOG[i]
  og=N0$OG[i]

  pathIn=paste(pathCDS, hog,"_",og,".aln.best.fas",sep="")
  
  if(file.exists(pathIn)){
    id.pauxi=N0[i,"Pauxi_pauxi"]

    if(id.pauxi%in%rownames(names)){
      name.pauxi=names[id.pauxi, "GeneName.Pauxi_pauxi"]

      if(length(grep("/",name.pauxi))!=0){
        name.pauxi=str_replace(name.pauxi, '/', '_')
      }
      
      if(!is.na(name.pauxi)){
        system(paste("cp ",pathIn, " ",pathOutput, "/",id.pauxi,"_",name.pauxi,".fa",sep=""))
      } else{
        system(paste("cp ",pathIn, " ",pathOutput, "/",id.pauxi,".fa",sep=""))
      }
    } else{
      system(paste("cp ",pathIn, " ",pathOutput, "/",id.pauxi,".fa",sep=""))
    }
  }

  if(i%%1000==0){
    print(i)
  }
}

###########################################################################



