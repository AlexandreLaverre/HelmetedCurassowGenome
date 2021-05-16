################################################################################

pathAnnot="../../results/genome_annotation/"
pathFigures="../../results/figures/"

assembly="MEGAHIT_RAGOUT"
braker="BRAKER_Ensembl103_multithread"
gemoma="GeMoMa/combined"

################################################################################

library(seqinr)
library(ade4)

################################################################################

## read protein sequences for all species, BRAKER2 data

files=system(paste("ls ",pathAnnot, assembly, "/", braker,"/OrthoFinder/ | grep fa$", sep=""), intern=T)
species=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))

################################################################################

proteins=list()
freqaa=list()

for(sp in species){
  print(sp)
  
  if(sp=="Pauxi_pauxi"){
    proteins[["BRAKER2"]]=read.fasta(paste(pathAnnot, assembly, "/", braker, "/OrthoFinder/", sp, ".fa", sep=""), seqtype="AA")
    freqaa[["BRAKER2"]]=as.numeric(table(factor(unlist(proteins[["BRAKER2"]]), levels=a())))
  } else{
    proteins[[sp]]=read.fasta(paste(pathAnnot, assembly, "/", braker, "/OrthoFinder/", sp, ".fa", sep=""), seqtype="AA")
    freqaa[[sp]]=as.numeric(table(factor(unlist(proteins[[sp]]), levels=a())))
  }  
}

################################################################################

## GeMoMa results

proteins[["GeMoMa"]]=read.fasta(paste(pathAnnot, assembly, "/",gemoma, "/final_predictions.faa", sep=""), seqtype="AA")
freqaa[["GeMoMa"]]=as.numeric(table(factor(unlist(proteins[["GeMoMa"]]), levels=a())))

################################################################################

freqaa=t(as.data.frame(freqaa))
colnames(freqaa)=a()
rownames(freqaa)=names(proteins)

################################################################################

afc <- dudi.coa(freqaa, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_BRAKER2_GeMoMa.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc$li[,1],afc$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc$eig[1]/sum(afc$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc$eig[2]/sum(afc$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc$li[,1]))/100
smally=diff(range(afc$li[,2]))/100

points(afc$li["BRAKER2",1],afc$li["BRAKER2",2], pch=20, col="red")
points(afc$li["GeMoMa",1],afc$li["GeMoMa",2], pch=20, col="blue")

legend("topleft", legend=c("BRAKER2", "GeMoMa"), pch=20, col=c("red", "blue"), inset=0.01)

dev.off()

################################################################################

## repeats for each annotation

repeats.braker2=read.table(paste(pathAnnot, assembly, "/", braker, "/overlap_repeats.txt", sep=""), h=T, stringsAsFactors=F)
repeats.gemoma=read.table(paste(pathAnnot, assembly, "/", gemoma, "/overlap_repeats.txt", sep=""), h=T, stringsAsFactors=F)

################################################################################
## select proteins with > 10% repeats, gemoma
repeats.gemoma=repeats.gemoma[which(repeats.gemoma$TranscriptID%in%names(proteins[["GeMoMa"]])),]
repeats.gemoma$Length=repeats.gemoma$End-repeats.gemoma$Start+1
txreplen.gemoma=tapply(repeats.gemoma$Length, as.factor(repeats.gemoma$TranscriptID), sum)
totlen.gemoma=unlist(lapply(proteins[["GeMoMa"]], length))*3 ## nucleotides
txreplen.gemoma.order=rep(0, length(totlen.gemoma))
names(txreplen.gemoma.order)=names(totlen.gemoma)
txreplen.gemoma.order[names(txreplen.gemoma)]=txreplen.gemoma
frrep.gemoma=txreplen.gemoma.order/totlen.gemoma

forbidden.gemoma=names(frrep.gemoma)[which(frrep.gemoma>0.1)]

## select proteins with > 10% repeats, braker
repeats.braker2=repeats.braker2[which(repeats.braker2$TranscriptID%in%names(proteins[["BRAKER2"]])),]
repeats.braker2$Length=repeats.braker2$End-repeats.braker2$Start+1
txreplen.braker2=tapply(repeats.braker2$Length, as.factor(repeats.braker2$TranscriptID), sum)
totlen.braker2=unlist(lapply(proteins[["BRAKER2"]], length))*3 ## nucleotides
txreplen.braker2.order=rep(0, length(totlen.braker2))
names(txreplen.braker2.order)=names(totlen.braker2)
txreplen.braker2.order[names(txreplen.braker2)]=txreplen.braker2
frrep.braker2=txreplen.braker2.order/totlen.braker2

forbidden.braker2=names(frrep.braker2)[which(frrep.braker2>0.1)]

################################################################################

proteins[["GEMOMA_norep"]]=proteins[["GeMoMa"]][setdiff(names(proteins[["GeMoMa"]]), forbidden.gemoma)]
proteins[["BRAKER2_norep"]]=proteins[["BRAKER2"]][setdiff(names(proteins[["BRAKER2"]]), forbidden.braker2)]

freqaa=rbind(freqaa, as.numeric(table(factor(unlist(proteins[["GeMoMa_norep"]]), levels=a()))))
rownames(freqaa)[nrow(freqaa)]="GeMoMa_norep"

freqaa=rbind(freqaa,as.numeric(table(factor(unlist(proteins[["BRAKER2_norep"]]), levels=a()))))
rownames(freqaa)[nrow(freqaa)]="BRAKER2_norep"

################################################################################

afc <- dudi.coa(freqaa, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_BRAKER2_GeMoMa_norepeats.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc$li[,1],afc$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc$eig[1]/sum(afc$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc$eig[2]/sum(afc$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc$li[,1]))/100
smally=diff(range(afc$li[,2]))/100

points(afc$li["BRAKER2",1],afc$li["BRAKER2",2], pch=20,col="red")
points(afc$li["GeMoMa",1],afc$li["GeMoMa",2], pch=20, col="blue")


points(afc$li["BRAKER2_norep",1],afc$li["BRAKER2_norep",2], pch=20, col="pink")
points(afc$li["GeMoMa_norep",1],afc$li["GeMoMa_norep",2], pch=20, col="steelblue")

legend("topleft", legend=c("BRAKER2", "GeMoMa", "BRAKER2, no repeats", "GeMoMa, no repeats"), pch=20, col=c("red", "blue", "pink", "steelblue"), inset=0.01)

dev.off()

################################################################################
