#########################################################################

pathResults="../../results/coding_gene_evolution/"

#########################################################################

library(ape)

#########################################################################

write.tree.adapted <- function (phy, file = "", append = FALSE, digits = 10, tree.names = FALSE) 
{
    if (!(inherits(phy, c("phylo", "multiPhylo")))) 
        stop("object \"phy\" has no trees")
    if (inherits(phy, "phylo")) 
        phy <- c(phy)
    N <- length(phy)
    res <- character(N)
    if (is.logical(tree.names)) {
        if (tree.names) {
            tree.names <- if (is.null(names(phy))) 
                character(N)
            else names(phy)
        }
        else tree.names <- character(N)
    }
    check_tips <- TRUE
    if (inherits(phy, "multiPhylo")) {
        if (!is.null(attr(phy, "TipLabel"))) {
            # attr(phy, "TipLabel") <- checkLabel(attr(phy, "TipLabel"))
            check_tips <- FALSE
        }
    }
    ## phy <- ape:::.uncompressTipLabel(phy)
    class(phy) <- NULL
    for (i in 1:N) res[i] <- ape:::.write.tree2(phy[[i]], digits = digits, 
        tree.prefix = tree.names[i], FALSE)
    if (file == "") 
        return(res)
    else cat(res, file = file, append = append, sep = "\n")
}

#########################################################################

helmeted=c("Anseranas_semipalmata", "Numida_meleagris", "Casuarius_casuarius", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

helmeted10=substr(helmeted, 1, 10)

rhinoceros=c("Bucorvus_abyssinicus", "Buceros_rhinoceros")
rhinoceros10=substr(rhinoceros,1,10)


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
  species=substr(species,1,10)

  this.tree=keep.tip(full.tree, species)
  this.tree$node.label=rep("", this.tree$Nnode)
  
  ## check if we need to label internal branch
  
  if(all(rhinoceros10%in%this.tree$tip.label)){
    anc=getMRCA(this.tree, tip=rhinoceros10)
    node.nb=anc-length(this.tree$tip.label)
    this.tree$node.label[node.nb]=" #1"

    stop()
  }

  ## add labels for external branches, helmeted birds
  
  this.helmeted=which(this.tree$tip.label%in%helmeted10)
  this.tree$tip.label[this.helmeted]=paste(this.tree$tip.label[this.helmeted], "#1", sep=" ")

  ## add labels for internal branches

  write.tree.adapted(this.tree,file=paste(pathResults, "CDS/",prefix,".branchmodel.tree",sep=""))

  nbdone=nbdone+1
  
  if(nbdone%%1000==0){
    print(paste("done ",nbdone,"trees"))
  }
}

#########################################################################

