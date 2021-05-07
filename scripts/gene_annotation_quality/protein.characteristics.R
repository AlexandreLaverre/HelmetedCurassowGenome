################################################################################

pathAnnot="../../results/genome_annotation/"
pathFigures="../../results/figures/"

assembly="MEGAHIT_RAGOUT"
method="BRAKER_Ensembl103_multithread"

################################################################################

library(seqinr)
library(ade4)

################################################################################

## read protein sequences for all species

files=system(paste("ls ",pathAnnot, assembly, "/", method, "/OrthoFinder/ | grep fa$", sep=""), intern=T)
species=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))

################################################################################

proteins=list()
propwithstop=list()
freqaa=list()
protlen=list()

for(sp in species){
  print(sp)
  
  proteins[[sp]]=read.fasta(paste(pathAnnot, assembly, "/", method, "/OrthoFinder/", sp, ".fa", sep=""), seqtype="AA")

  nbstop=unlist(lapply(proteins[[sp]], function(x) length(which(x=="*"))))
  len=unlist(lapply(proteins[[sp]], length))

  protlen[[sp]]=len
  propwithstop[[sp]]=100*length(which(nbstop>0))/length(proteins[[sp]])
  freqaa[[sp]]=as.numeric(table(factor(unlist(proteins[[sp]]), levels=a())))
}

freqaa=t(as.data.frame(freqaa))
colnames(freqaa)=a()
rownames(freqaa)=species

################################################################################

afc <- dudi.coa(freqaa, scann = FALSE, nf = 5)

################################################################################

## premier plan factoriel

pdf(paste(pathFigures, "COA_AminoAcidComposition_BRAKER2.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc$li[,1],afc$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc$eig[1]/sum(afc$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc$eig[2]/sum(afc$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc$li[,1]))/100
smally=diff(range(afc$li[,2]))/100

text(labels=rownames(afc$li)[50], x=afc$li[50,1]+smallx, y=afc$li[50,2]+smally, adj=0)

dev.off()

################################################################################

## length distribution

pdf(paste(pathFigures, "ProteinLength_BRAKER2.pdf",sep=""), width=4, height=11)
par(mar=c(1,3.1, 2.1, 1.1))

boxplot(protlen, horizontal=T)

dev.off()

################################################################################
