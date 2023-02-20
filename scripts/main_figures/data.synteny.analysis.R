######################################################################

pathResults="../../results/synteny_analysis/"
pathFigures="../../results/figures/"

######################################################################

species=c("Basiliscus_vittatus", "Pauxi_pauxi")
otherspecies=c("Gallus_gallus", "Anolis_carolinensis")

######################################################################

synteny=list()

for(i in 1:2){
    sp=species[i]
    othersp=otherspecies[i]

    synteny[[sp]]=read.table(paste(pathResults, sp, "/OrthoGeneCoordinates_",othersp,".txt",sep=""),h=T, stringsAsFactors=F)
}

######################################################################

save(synteny, file=paste(pathFigures, "RData/synteny.RData",sep=""))

######################################################################
