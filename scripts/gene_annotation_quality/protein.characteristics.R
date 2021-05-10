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

## COA, all proteins

afc <- dudi.coa(freqaa, scann = FALSE, nf = 5)

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

par(mar=c(4,8.1, 2.1, 1.1))
par(mgp=c(3, 0.5, 0))
par(cex.axis=0.85)
boxplot(protlen, horizontal=T, names=rep("", length(protlen)), outline=F, xlim=c(0.5, length(protlen)+0.5), xaxs="i")
mtext(names(protlen), side=2, at=1:length(protlen), las=2, line=1, cex=0.5, font=3)

mtext("protein length", side=1, line=2, cex=0.95)

dev.off()

################################################################################

filtered.prot=list()
freqaa.filtered=list()
minlen=100

for(sp in species){
  print(sp)

  nbstop=unlist(lapply(proteins[[sp]], function(x) length(which(x=="*"))))
  len=protlen[[sp]]

  ok=which(len>=minlen & nbstop==0)
  filtered.prot[[sp]]=proteins[[sp]][ok]
  freqaa.filtered[[sp]]=as.numeric(table(factor(unlist(filtered.prot[[sp]]), levels=a())))
}

freqaa.filtered=t(as.data.frame(freqaa.filtered))
colnames(freqaa.filtered)=a()
rownames(freqaa.filtered)=species

################################################################################

## COA, filtered proteins

afc.filtered <- dudi.coa(freqaa.filtered, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_BRAKER2_MinLength100_NoStop.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc.filtered$li[,1],afc.filtered$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc.filtered$eig[1]/sum(afc.filtered$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc.filtered$eig[2]/sum(afc.filtered$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc.filtered$li[,1]))/100
smally=diff(range(afc.filtered$li[,2]))/100

text(labels=rownames(afc.filtered$li)[50], x=afc.filtered$li[50,1]+smallx, y=afc.filtered$li[50,2]+smally, adj=0)

dev.off()

################################################################################

orthogroups=read.table(paste(pathAnnot, assembly, "/", method, "/OrthoFinder/OrthoFinder/Results_May06_2/Orthogroups/Orthogroups.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

inortho=apply(orthogroups, 2, function(x) unlist(lapply(x, function(y) unlist(strsplit(y, split=", ")))))

## take proteins in orthogroups

ortho.prot=list()
freqaa.ortho=list()

for(sp in species){
  print(sp)

  ok=intersect(names(proteins[[sp]]), inortho[[sp]])
  ortho.prot[[sp]]=proteins[[sp]][ok]
  freqaa.ortho[[sp]]=as.numeric(table(factor(unlist(ortho.prot[[sp]]), levels=a())))
}

freqaa.ortho=t(as.data.frame(freqaa.ortho))
colnames(freqaa.ortho)=a()
rownames(freqaa.ortho)=species

################################################################################
## COA, ortho proteins

afc.ortho <- dudi.coa(freqaa.ortho, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_BRAKER2_InOrthogroups.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc.ortho$li[,1],afc.ortho$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc.ortho$eig[1]/sum(afc.ortho$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc.ortho$eig[2]/sum(afc.ortho$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc.ortho$li[,1]))/100
smally=diff(range(afc.ortho$li[,2]))/100

text(labels=rownames(afc.ortho$li)[50], x=afc.ortho$li[50,1]+smallx, y=afc.ortho$li[50,2]+smally, adj=0)

dev.off()

################################################################################

## multi-exonic genes

multiex=readLines(paste(pathAnnot, assembly, "/", method, "/multiexonic_genes_proteins.txt", sep=""))

all.filters=intersect(multiex, inortho[["Pauxi_pauxi"]]) ## multiexonic, orthogroups
all.filters=intersect(all.filters, names(filtered.prot[["Pauxi_pauxi"]])) ## min length, no stop

freqaa.allfilters=freqaa.ortho
freqaa.allfilters["Pauxi_pauxi",]=as.numeric(table(factor(unlist(proteins[["Pauxi_pauxi"]][all.filters]), levels=a())))

## COA, all filters

afc.allfilters <- dudi.coa(freqaa.allfilters, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_BRAKER2_MinLength100_NoStop_InOrthogroups_Multiexonic.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc.allfilters$li[,1],afc.allfilters$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc.allfilters$eig[1]/sum(afc.allfilters$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc.allfilters$eig[2]/sum(afc.allfilters$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc.allfilters$li[,1]))/100
smally=diff(range(afc.allfilters$li[,2]))/100

text(labels=rownames(afc.allfilters$li)[50], x=afc.allfilters$li[50,1]+smallx, y=afc.allfilters$li[50,2]+smally, adj=0)

dev.off()

################################################################################

propaa.allfilters=freqaa.allfilters/apply(freqaa.allfilters,1,sum)

################################################################################
