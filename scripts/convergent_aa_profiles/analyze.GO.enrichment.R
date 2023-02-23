#########################################################################################################

pathAnnot="../../data/ensembl_annotations/Human/"
pathPelican="../../results/coding_gene_evolution/all_species"

dataset="birds"

#########################################################################################################

cc=system(paste("grep cellular_component$ ", pathAnnot, "GeneOntology_Ensembl103.txt",sep=""), intern=T)

cc.gene=unlist(lapply(cc, function(x) unlist(strsplit(x, split="\t"))[1]))
cc.cat=unlist(lapply(cc, function(x) unlist(strsplit(x, split="\t"))[3]))

#########################################################################################################

all.sites=

#########################################################################################################
