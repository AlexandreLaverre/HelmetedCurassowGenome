#########################################################################

pathResults="../../results/coding_gene_evolution/"

#########################################################################

helmeted=c("Anseranas_semipalmata", "Numida_meleagris", "Casuarius_casuarius", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

helmeted10=substr(helmeted, 1, 10)

#########################################################################

library(ape)

#########################################################################

full.tree=read.tree(paste(pathResults, "species_tree_rooted.txt",sep=""))
full.tree$node.label <- NULL

full.tree$tip.label=substr(full.tree$tip.label,1,10)

#########################################################################

files=system(paste("ls ",pathResults, "CDS/ | grep unaln",sep=""), intern=T)

#########################################################################

nbdone=0

for(file in files){
  prefix=paste(unlist(strsplit(file, split="\\."))[1:2], collapse=".")

  species=system(paste("grep \">\" ",pathResults,"/CDS/",file,sep=""), intern=T)
  species=unlist(lapply(species, function(x) substr(x,2,nchar(x))))

  this.tree=keep.tip(full.tree, species)
  this.tree$node.label=rep("", this.tree$Nnode)
  
  ## check if we need to label internal branch
  
  if(all(c("Bucorvus_abyssinicus", "Buceros_rhinoceros")%in%this.tree$tip.label)){
    anc=getMRCA(this.tree, tip=c("Bucorvus_abyssinicus", "Buceros_rhinoceros"))
    node.nb=anc-length(this.tree$tip.label)
    this.tree$node.label[node.nb]="#1" 
  }

  ## add labels for external branches, helmeted birds
  
  this.helmeted=which(this.tree$tip.label%in%helmeted)
  this.tree$tip.label[this.helmeted]=paste(this.tree$tip.label[this.helmeted], "#1", sep=" ")

  ## add labels for internal branches

  write.tree(this.tree,file=paste(pathResults, "CDS/",prefix,".branchmodel.tree",sep=""))

  nbdone=nbdone+1
  
  if(nbdone%%1000==0){
    print(paste("done ",nbdone,"trees"))
  }
}

#########################################################################

