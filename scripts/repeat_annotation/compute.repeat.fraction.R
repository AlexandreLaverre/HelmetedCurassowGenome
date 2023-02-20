######################################################################

pathUCSC="../../data/RepeatMasker/"
pathInHouse="../../results/repeats/"

######################################################################

species=c("Chicken", "Lizard", "Basiliscus_vittatus", "Pauxi_pauxi")
annots=c("UCSC", "UCSC", "InHouse", "InHouse")
files=c("RepeatMasker_UCSC_galGal6.txt", "RepeatMasker_UCSC_anoCar2.txt", "CombinedRepeatAnnotations.txt", "CombinedRepeatAnnotations.txt")

######################################################################

for(i in 1:length(species)){
    sp=species[i]
    annot=annots[i]
    file=files[i]

    print(sp)

    if(annot=="UCSC"){
      repeats=read.table(paste(pathUCSC, sp, "/", file, sep=""), h=F, stringsAsFactors=F,sep="\t", quote="\"")
      colnames(repeats)[12]="Class"
      colnames(repeats)[6]="Chr"
      colnames(repeats)[7]="Start"
      colnames(repeats)[8]="End"
      repeats$Length=repeats$End-repeats$Start

      repeats$Class[which(repeats$Class=="DNA?")]="DNA"
      repeats$Class[which(repeats$Class=="RC?")]="RC"
      repeats$Class[which(repeats$Class=="LTR?")]="LTR"
      repeats$Class[which(repeats$Class=="LINE?")]="LINE"
      repeats$Class[which(repeats$Class=="SINE")]="SINE"

      pathOutput=paste(pathUCSC, sp, "/TotalLengthByClass.txt", sep="")
    }

    if(annot=="InHouse"){
      repeats=read.table(paste(pathInHouse, sp, "/MEGAHIT_RAGOUT/RepeatMasker/", file, sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
      repeats$Class=unlist(lapply(repeats$Class.Family, function(x) unlist(strsplit(x, split="\\/"))[1]))
      repeats$Length=repeats$End-repeats$Start+1

      repeats$Class[which(repeats$Class=="DNA?")]="DNA"
      repeats$Class[which(repeats$Class=="RC?")]="RC"
      repeats$Class[which(repeats$Class=="LTR?")]="LTR"
      repeats$Class[which(repeats$Class=="LINE?")]="LINE"
      repeats$Class[which(repeats$Class=="SINE")]="SINE"

      pathOutput=paste(pathInHouse, sp, "/MEGAHIT_RAGOUT/RepeatMasker/TotalLengthByClass.txt", sep="")
    }

    tot.len=tapply(repeats$Length, as.factor(repeats$Class), sum)

    df=data.frame("RepeatClass"=levels(as.factor(repeats$Class)),"TotalLength"=tot.len)

    write.table(df, file=pathOutput, row.names=F, col.names=T, sep="\t", quote=F)

}

######################################################################
