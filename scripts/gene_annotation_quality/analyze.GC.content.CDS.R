##############################################################################

library(seqinr)
library(ade4)

#############################################################################

pathEnsembl="../../data/coding_sequences/Ensembl103/primary_transcripts/"
pathAnnot="../../results/genome_annotation/"
pathFigures="../../results/figures/"

################################################################################

ensbirds=c("Accipiter_nisus",  "Amazona_collaria", "Anas_platyrhynchos", "Anas_platyrhynchos_platyrhynchos", "Anas_zonorhyncha","Anser_brachyrhynchus", "Anser_cygnoides", "Apteryx_haastii", "Apteryx_owenii", "Apteryx_rowi", "Aquila_chrysaetos_chrysaetos", "Athene_cunicularia", "Bubo_bubo", "Buteo_japonicus", "Cairina_moschata_domestica", "Calidris_pugnax", "Calidris_pygmaea", "Camarhynchus_parvulus", "Catharus_ustulatus", "Corvus_moneduloides", "Coturnix_japonica", "Chrysolophus_pictus",  "Cyanistes_caeruleus", "Dromaius_novaehollandiae", "Erythrura_gouldiae", "Falco_tinnunculus", "Ficedula_albicollis", "Gallus_gallus", "Geospiza_fortis",  "Junco_hyemalis", "Lepidothrix_coronata", "Lonchura_striata_domestica", "Malurus_cyaneus_samueli", "Manacus_vitellinus", "Meleagris_gallopavo", "Melopsittacus_undulatus",  "Nothoprocta_perdicaria", "Numida_meleagris", "Otus_sunia", "Parus_major",  "Pavo_cristatus",  "Phasianus_colchicus", "Serinus_canaria",  "Stachyris_ruficeps", "Strigops_habroptila", "Strix_occidentalis_caurina", "Struthio_camelus_australis", "Taeniopygia_guttata",   "Zonotrichia_albicollis", "Zosterops_lateralis_melanops")

mammals=c("Mus_musculus", "Homo_sapiens")

squamates=c("Anolis_carolinensis","Varanus_komodoensis", "Naja_naja", "Notechis_scutatus","Pseudonaja_textilis", "Salvator_merianae", "Podarcis_muralis", "Pogona_vitticeps","Laticauda_laticaudata")

crocodile=c("Crocodylus_porosus")

turtles=c("Terrapene_carolina_triunguis","Chelonoidis_abingdonii", "Chelydra_serpentina", "Chrysemys_picta_bellii", "Gopherus_agassizii", "Gopherus_evgoodei", "Pelodiscus_sinensis", "Pelusios_castaneus")

tuatara=c("Sphenodon_punctatus")

ensemblspecies=c(ensbirds, mammals, squamates, crocodile, turtles, tuatara)
ensemblpaths=system(paste("ls ",pathEnsembl, "| grep -v dmnd"), intern=T)

################################################################################

## new species

newspecies=system(paste("ls ",pathAnnot, sep=""), intern=T)
assemblies=unlist(lapply(newspecies, function(x) system(paste("ls ",pathAnnot,"/",x, sep=""), intern=T)))
names(assemblies)=newspecies

otherbirds=setdiff(newspecies, c("Basiliscus_vittatus", "Pauxi_pauxi"))
allbirds=c(ensbirds, otherbirds)

species=c(ensemblspecies, newspecies)

################################################################################

cds=list()

for(sp in ensbirds){
  print(sp)

  path=grep(paste(sp, "\\.",sep=""), ensemblpaths, value=T)
  cds[[sp]]=read.fasta(paste(pathEnsembl, path, sep=""), seqtype="DNA", forceDNAtolower=FALSE)

  cds.gc3=unlist(lapply(cds[[sp]], GC3))
  cds.gc=unlist(lapply(cds[[sp]], GC))

  xlim=c(0,1)
  d.gc3=density(cds.gc3, bw=0.015)
  d.gc=density(cds.gc3, bw=0.015)
  ylim=range(c(d.gc3$y, d.gc$y))

  pdf(file=paste(pathFigures, "GCContent_CDS_",sp,"_Ensembl103.pdf",sep=""), width=6,height=4)
  plot(d.gc3, xlim=xlim, ylim=ylim, col="red", xlab="GC content", ylab="density", main=sp)
  lines(d.gc, col="black")
  legend("topright",lty=1, col=c("red", "black"), legend=c("GC3", "GC"), inset=0.01)
  dev.off()
}

for(parset in c("filtered_predictions", "filtered_predictions_minDiamondProteinFraction0.25_minLength100_maxFractionRepeats0.25", "filtered_predictions_minDiamondProteinFraction0.25_minLength100_maxFractionRepeats0.5", "filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.25","filtered_predictions_minDiamondProteinFraction0.25_minLength70_maxFractionRepeats0.5")){
    print(parset)

    for(sp in otherbirds){
        print(sp)

        assembly=assemblies[sp]
        cds[[sp]]=read.fasta(paste(pathAnnot, sp, "/", assembly, "/GeMoMa/combined/",parset,".cds.fa", sep=""), seqtype="DNA", forceDNAtolower=FALSE)

        cds.gc3=unlist(lapply(cds[[sp]], GC3))
        cds.gc=unlist(lapply(cds[[sp]], GC))

        xlim=c(0,1)
        d.gc3=density(cds.gc3, bw=0.015)
        d.gc=density(cds.gc3, bw=0.015)
        ylim=range(c(d.gc3$y, d.gc$y))

        pdf(file=paste(pathFigures, "GCContent_CDS_",sp,"_",parset,".pdf",sep=""), width=6,height=4)
        plot(d.gc3, xlim=xlim, ylim=ylim, col="red", xlab="GC content", ylab="density", main=sp)
        lines(d.gc, col="black")
        legend("topright",lty=1, col=c("red", "black"), legend=c("GC3", "GC"), inset=0.01)
        dev.off()
    }
}

#############################################################################
