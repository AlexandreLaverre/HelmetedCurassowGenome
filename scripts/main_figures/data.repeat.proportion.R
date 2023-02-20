######################################################################

pathUCSC="../../data/RepeatMasker/"
pathInHouse="../../results/repeats/"
pathFigures="../../results/figures/"

######################################################################

species=c("Chicken", "Lizard", "Basiliscus_vittatus", "Pauxi_pauxi")
annots=c("UCSC", "UCSC", "InHouse", "InHouse")
synspecies=c("Gallus_gallus", "Anolis_carolinensis", "Basiliscus_vittatus", "Pauxi_pauxi")

######################################################################

repfr=list()

for(i in 1:length(species)){
    sp=species[i]
    annot=annots[i]

    if(annot=="UCSC"){
        pathInput=paste(pathUCSC, sp, "/TotalLengthByClass.txt", sep="")
    }

    if(annot=="InHouse"){
        pathInput=paste(pathInHouse, sp, "/MEGAHIT_RAGOUT/RepeatMasker/TotalLengthByClass.txt", sep="")
    }

    repfr[[sp]]=read.table(pathInput, h=T, stringsAsFactors=F)

}

######################################################################

save(repfr, file=paste(pathFigures, "RData/repeat.fraction.RData",sep=""))

######################################################################
