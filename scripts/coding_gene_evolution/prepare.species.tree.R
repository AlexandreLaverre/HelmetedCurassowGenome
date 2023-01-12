########################################################################

library(ape)

########################################################################

pathResults="../../results/coding_gene_evolution/"

########################################################################

for(spset in c("all_species", "without_chameleons")){

  if(spset=="without_chameleons"){
    pathOrthoFinder="../../results/gene_families/OrthoFinder/without_chameleons/iqtree/Results_Sep21/"
  }

  if(spset=="all_species"){
    pathOrthoFinder="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/"
  }

  tree=read.tree(paste(pathOrthoFinder, "/Species_Tree/SpeciesTree_rooted.txt",sep=""))
  
  for(dataset in c("all_species", "birds", "squamates")){
    
    if(spset=="without_chameleons" & dataset=="all_species"){
      rooted=root(tree, outgroup=c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis"))
    }

    if(spset=="all_species" & dataset=="all_species"){
      rooted=root(tree, outgroup=c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis", "Chamaeleo_chamaeleon_recticrista","Chamaeleo_calyptratus"))
    }
    
    if(dataset=="birds"){
      filtered=drop.tip(tree, intersect(c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis", "Chamaeleo_chamaeleon_recticrista","Chamaeleo_calyptratus"), tree$tip.label))
      rooted=root(filtered, outgroup=c("Struthio_camelus_australis", "Casuarius_casuarius", "Dromaius_novaehollandiae"))
    }

    if(dataset=="squamates"){
      filtered=keep.tip(tree, intersect(c("Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps", "Basiliscus_vittatus", "Anolis_carolinensis", "Naja_naja", "Notechis_scutatus", "Pseudonaja_textilis", "Chamaeleo_chamaeleon_recticrista","Chamaeleo_calyptratus"), tree$tip.label))
      rooted=root(filtered, outgroup=c("Pseudonaja_textilis", "Notechis_scutatus", "Naja_naja"))
    }
    
  
  }

  write.tree(rooted, file=paste(pathResults,spset, "/", dataset, "/species_tree.txt",sep=""))
}

########################################################################
