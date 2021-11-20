#########################################################################

pathResults="../../results/coding_gene_evolution/"

#########################################################################

library(ape)

#########################################################################

full.tree=read.tree(paste(pathResults, "species_tree_rooted.txt",sep=""))

#########################################################################

files=system(paste("ls ",pathResults, "CDS/",sep=""), intern=T)

#########################################################################

for(file in files){
  prefix=paste(unlist(strsplit(file, split="\\."))[1:2], collapse=".")

  species=system(paste("grep \">\" ",pathResults,"/CDS/",file,sep=""), intern=T)
  species=unlist(lapply(species, function(x) substr(x,2,nchar(x))))

  this.tree=keep.tip(full.tree, species)

  write.tree(this.tree,file=paste(pathResults, "CDS/",prefix,".tree",sep=""))
}

#########################################################################

