#########################################################################

pathResults="../../results/coding_gene_evolution/"

#########################################################################

library(ape)

#########################################################################

for(dataset in c("all_species", "birds")){

    full.tree=read.tree(paste(pathResults, dataset, "/species_tree.txt",sep=""))
    full.tree$node.label <- NULL

    #########################################################################

    files=system(paste("ls ",pathResults, dataset, "/CDS/ | grep unaln.fa",sep=""), intern=T)

    #########################################################################

    for(file in files){
        prefix=paste(unlist(strsplit(file, split="\\."))[1:2], collapse=".")

        species=system(paste("grep \">\" ",pathResults, dataset, "/CDS/",file,sep=""), intern=T)
        species=unlist(lapply(species, function(x) substr(x,2,nchar(x))))

        this.tree=keep.tip(full.tree, species)

        write.tree(this.tree,file=paste(pathResults, dataset, "/CDS/",prefix,".tree",sep=""))
    }
}

#########################################################################
