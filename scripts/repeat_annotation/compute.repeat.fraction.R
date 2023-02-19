######################################################################

pathUCSC="../../data/RepeatMasker/"
pathInHouse="../../results/repeats/"

######################################################################

species=c("Chicken", "Lizard", "Basiliscus_vittatus", "Pauxi_pauxi")
sources=c("UCSC", "UCSC", "InHouse", "InHouse")
files=c("RepeatMasker_UCSC_galGal6.txt", "RepeatMasker_UCSC_anoCar2.txt", "CombinedRepeatAnnotations.txt", "CombinedRepeatAnnotations.txt")

######################################################################

for(i in 1:length(species)){
    sp=species[i]
    source=sources[i]
    file=files[i]

    print(sp)

    if(source=="UCSC"){
        repeats=read.table(paste(pathUCSC, sp, "/", file, sep=""), h=T, stringsAsFactors=F,sep="\t", quote="\"")
    }

    if(source=="InHouse"){
        repeats=read.table(paste(pathInHouse, sp, "/MEGAHIT_RAGOUT/RepeatMasker/", file, sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    }
}

######################################################################
