#########################################################################

pathResults="../../results/coding_gene_evolution/"

#########################################################################

library(ape)

#########################################################################

for(spset in c("all_species", "without_chameleons")){
  for(dataset in c("all_species", "birds")){
    
    full.tree=read.tree(paste(pathResults,spset,"/", dataset, "/species_tree.txt",sep=""))
    full.tree$node.label <- NULL
    
#########################################################################
    
    files=system(paste("ls ",pathResults, spset, "/", dataset, "/CDS/ | grep unaln.fa",sep=""), intern=T)
    
#########################################################################
    
    for(file in files){
      prefix=paste(unlist(strsplit(file, split="\\."))[1:2], collapse=".")
      
      species=system(paste("grep \">\" ",pathResults, spset, "/", dataset, "/CDS/",file,sep=""), intern=T)
      species=unlist(lapply(species, function(x) substr(x,2,nchar(x))))
      
      this.tree=keep.tip(full.tree, species)
      
      write.tree(this.tree,file=paste(pathResults, spset, "/", dataset, "/CDS/",prefix,".tree",sep=""))
    }
}

#########################################################################
