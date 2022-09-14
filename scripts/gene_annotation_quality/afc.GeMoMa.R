################################################################################

pathEnsembl="../../data/protein_sequences/Ensembl103/primary_transcripts/"
pathAnnot="../../results/genome_annotation/"
pathFigures="../../results/figures/"

################################################################################

library(seqinr)
library(ade4)

################################################################################

ensbirds=c("Accipiter_nisus", "Alectura_lathami", "Amazona_collaria", "Anas_platyrhynchos", "Anas_platyrhynchos_platyrhynchos", "Anas_zonorhyncha","Anser_brachyrhynchus", "Anser_cygnoides", "Apteryx_haastii", "Apteryx_owenii", "Apteryx_rowi", "Aquila_chrysaetos_chrysaetos", "Athene_cunicularia", "Bubo_bubo", "Buteo_japonicus", "Cairina_moschata_domestica", "Calidris_pugnax", "Calidris_pygmaea", "Camarhynchus_parvulus", "Casuarius_casuarius", "Catharus_ustulatus", "Corvus_moneduloides", "Coturnix_japonica", "Chrysolophus_pictus",  "Cyanistes_caeruleus", "Dromaius_novaehollandiae", "Erythrura_gouldiae", "Falco_tinnunculus", "Ficedula_albicollis", "Gallus_gallus", "Geospiza_fortis",  "Junco_hyemalis", "Lepidothrix_coronata", "Lonchura_striata_domestica", "Malurus_cyaneus_samueli", "Manacus_vitellinus", "Meleagris_gallopavo", "Melopsittacus_undulatus",  "Nothoprocta_perdicaria", "Numida_meleagris", "Otus_sunia", "Parus_major",  "Pavo_cristatus", "Penelope_pileata", "Phasianus_colchicus", "Serinus_canaria",  "Stachyris_ruficeps", "Strigops_habroptila", "Strix_occidentalis_caurina", "Struthio_camelus_australis", "Taeniopygia_guttata",   "Zonotrichia_albicollis", "Zosterops_lateralis_melanops")

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

################################################################################

species=c(ensemblspecies, newspecies)

################################################################################

proteins=list()

for(sp in species){
  print(sp)

  if(sp%in%newspecies){
      assembly=assemblies[sp]
      proteins[[sp]]=read.fasta(paste(pathAnnot, sp, "/", assembly, "/GeMoMa/combined/final_annotations.faa", sep=""), seqtype="AA")

  } else{
      path=grep(paste(sp, "\\.",sep=""), ensemblpaths, value=T)
      proteins[[sp]]=read.fasta(paste(pathEnsembl, path, sep=""), seqtype="AA")

  }
}

################################################################################

freqaa=list()

for(sp in species){
    print(sp)
    freqaa[[sp]]=as.numeric(table(factor(unlist(proteins[[sp]]), levels=a())))
}

################################################################################

freqaa=t(as.data.frame(freqaa))
colnames(freqaa)=a()
rownames(freqaa)=names(proteins)

freqaa=freqaa[,which(colnames(freqaa)!="*")]

################################################################################

print("AFC")

afc <- dudi.coa(freqaa, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_taxonomy.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc$li[,1],afc$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1, type="n")

mtext(paste0("F1 (", round(100*afc$eig[1]/sum(afc$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc$eig[2]/sum(afc$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc$li[,1]))/100
smally=diff(range(afc$li[,2]))/100

points(afc$li[ensbirds,1],afc$li[ensbirds, 2], pch=21, col="indianred", bg="indianred", cex=1.2)
points(afc$li[otherbirds,1],afc$li[otherbirds, 2], pch=21, col="brown", bg="brown", cex=1.2)
points(afc$li[crocodile,1],afc$li[crocodile,2], pch=21, col="gray40", bg="gray40", cex=1.2)
points(afc$li[turtles,1],afc$li[turtles,2], pch=21, col="slateblue", bg="steelblue", cex=1.2)
points(afc$li[squamates,1],afc$li[squamates,2], pch=21, col="seagreen", bg="seagreen", cex=1.2)
points(afc$li[tuatara,1],afc$li[tuatara,2], pch=21, col="green", bg="green", cex=1.2)
points(afc$li[mammals,1],afc$li[mammals,2], pch=21, col="orange", bg="orange", cex=1.2)

points(afc$li["Pauxi_pauxi",1],afc$li["Pauxi_pauxi",2], pch=21, bg="red", col="red", cex=1.1)
points(afc$li["Basiliscus_vittatus",1],afc$li["Basiliscus_vittatus",2], pch=21, bg="blue", col="blue", cex=1.1)

legend("topleft", legend=c("helmeted curassow", "brown basilisc", "Ensembl birds", "GeMoMa birds", "crocodile", "turtles", "squamates",  "tuatara", "mammals"), pch=21, pt.bg=c("red", "blue", "indianred", "brown", "gray40", "steelblue", "seagreen", "green", "orange"), col=c("red", "blue", "indianred", "brown", "gray40", "steelblue", "seagreen", "green", "orange"), inset=0.01)

dev.off()

################################################################################

