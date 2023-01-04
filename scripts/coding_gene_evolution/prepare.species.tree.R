########################################################################

library(ape)

########################################################################

pathOrthoFinder="../../results/gene_families/OrthoFinder/iqtree/Results_Sep21/"
pathResults="../../results/coding_gene_evolution/"

########################################################################

tree=read.tree(paste(pathOrthoFinder, "/Species_Tree/SpeciesTree_rooted.txt",sep=""))

########################################################################

for(dataset in c("all_species", "birds")){

    if(dataset=="all_species"){
        rooted=root(tree, outgroup=c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis"))
    }

    if(dataset=="birds"){
        filtered=drop.tip(tree, c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis"))
        rooted=root(filtered, outgroup=c("Struthio_camelus_australis", "Casuarius_casuarius", "Dromaius_novaehollandiae"))
    }

    write.tree(rooted, file=paste(pathResults,dataset, "/species_tree.txt",sep=""))
}

########################################################################
