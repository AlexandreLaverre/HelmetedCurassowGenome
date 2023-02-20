######################################################################

pathFigures="../../results/figures/"

######################################################################

gcinfo=list()

species=c("Basiliscus_vittatus", "Pauxi_pauxi")

for(i in 1:2){
    sp=species[i]

    pathGenome=paste("../../results/genome_assembly/", sp, "/MEGAHIT_RAGOUT/",sep="")

    gcsize=read.table(paste(pathGenome, "GCContent_ChromosomeSize.txt",sep=""), h=T, stringsAsFactors=F)

    gc3=read.table(paste(pathGenome, "MeanGC3Content.txt",sep=""), h=T, stringsAsFactors=F)
    rownames(gc3)=gc3$Chr

    gcsize$MeanGC3=rep(NA, nrow(gcsize))
    gcsize[which(gcsize$Chr%in%gc3$Chr), "MeanGC3"]=gc3[gcsize$Chr[which(gcsize$Chr%in%gc3$Chr)], "MeanGC3"]

    gcinfo[[sp]]=gcsize
}

######################################################################

save(gcinfo, file=paste(pathFigures, "RData/gccontent.RData",sep=""))

######################################################################
