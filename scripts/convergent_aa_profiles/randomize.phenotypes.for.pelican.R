###########################################################################

library(ape)

###########################################################################

all=c("Struthio_camelus_australis", "Casuarius_casuarius", "Dromaius_novaehollandiae", "Gallus_gallus", "Meleagris_gallopavo", "Numida_meleagris", "Penelope_pileata", "Alectura_lathami", "Anseranas_semipalmata", "Anas_platyrhynchos_platyrhynchos", "Anser_cygnoides", "Balearica_regulorum", "Grus_americana", "Calidris_pugnax", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Upupa_epops", "Rhinopomastus_cyanomelas", "Strix_occidentalis_caurina", "Ficedula_albicollis", "Anser_brachyrhynchus", "Aquila_chrysaetos_chrysaetos", "Geospiza_fortis", "Parus_major", "Serinus_canaria", "Strigops_habroptila", "Coturnix_japonica", "Anolis_carolinensis", "Podarcis_muralis", "Pogona_vitticeps", "Pseudonaja_textilis", "Salvator_merianae", "Naja_naja", "Notechis_scutatus", "Chamaeleo_chamaeleon_recticrista", "Chamaeleo_calyptratus", "Basiliscus_vittatus", "Pauxi_pauxi")

withoutchameleons=c("Struthio_camelus_australis", "Casuarius_casuarius", "Dromaius_novaehollandiae", "Gallus_gallus", "Meleagris_gallopavo", "Numida_meleagris", "Penelope_pileata", "Alectura_lathami", "Anseranas_semipalmata", "Anas_platyrhynchos_platyrhynchos", "Anser_cygnoides", "Balearica_regulorum", "Grus_americana", "Calidris_pugnax", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Upupa_epops", "Rhinopomastus_cyanomelas", "Strix_occidentalis_caurina", "Ficedula_albicollis", "Anser_brachyrhynchus", "Aquila_chrysaetos_chrysaetos", "Geospiza_fortis", "Parus_major", "Serinus_canaria", "Strigops_habroptila", "Coturnix_japonica", "Anolis_carolinensis", "Podarcis_muralis", "Pogona_vitticeps", "Pseudonaja_textilis", "Salvator_merianae", "Naja_naja", "Notechis_scutatus",  "Basiliscus_vittatus", "Pauxi_pauxi")

birds=c("Struthio_camelus_australis", "Casuarius_casuarius", "Dromaius_novaehollandiae", "Gallus_gallus", "Meleagris_gallopavo", "Numida_meleagris", "Penelope_pileata", "Alectura_lathami", "Anseranas_semipalmata", "Anas_platyrhynchos_platyrhynchos", "Anser_cygnoides", "Balearica_regulorum", "Grus_americana", "Calidris_pugnax", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Upupa_epops", "Rhinopomastus_cyanomelas", "Strix_occidentalis_caurina", "Ficedula_albicollis", "Anser_brachyrhynchus", "Aquila_chrysaetos_chrysaetos", "Geospiza_fortis", "Parus_major", "Serinus_canaria", "Strigops_habroptila", "Coturnix_japonica",  "Pauxi_pauxi")

###########################################################################

set.seed(19)

protuberance=c(sample(birds, size=8, replace=F), sample(setdiff(all, birds), size=2, replace=F))

###########################################################################

for(spset in c("all_species")){ ## other option is "without_chameleons"
    for(dataset in c("all_species", "birds", "squamates")){

        pathResults=paste("../../results/coding_gene_evolution/",spset, "/",dataset, sep="")

        full.tree=read.tree(paste(pathResults, "/species_tree.txt",sep=""))

        if(file.exists(paste(pathResults, "/species_tree_nobootstrap.txt",sep=""))){
            print("tree already there")
        }  else {
            full.tree$node.label <- NULL
            write.tree(full.tree, paste(pathResults, "/species_tree_nobootstrap.txt",sep=""))
        }

        for(phenotype in c("random")){

            all.species=full.tree$tip.label

            phenotype.data=data.frame("species"=all.species, "trait"=rep(NA, length(all.species)))
            phenotype.data$trait[which(!(phenotype.data$species%in%protuberance))]="no_protuberance"

            if(phenotype == "random"){
                phenotype.data$trait[which(phenotype.data$species%in%protuberance)]="protuberance"
            } else{
                splist=get(phenotype, envir=globalenv())
                phenotype.data$trait[which(phenotype.data$species%in%splist)]=phenotype
                phenotype.data$trait[which(phenotype.data$species%in%setdiff(protuberance, splist))]="other_protuberance"
            }

            write.table(phenotype.data, file=paste(pathResults, "/phenotype_data_",phenotype,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
        }
    }
}

###########################################################################
