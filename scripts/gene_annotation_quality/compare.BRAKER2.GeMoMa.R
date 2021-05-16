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
  
  proteins[[sp]]=read.fasta(paste(pathAnnot, assembly, "/", braker, "/OrthoFinder/", sp, ".fa", sep=""), seqtype="AA")

  
  if(sp=="Pauxi_pauxi"){
    freqaa[["BRAKER2"]]=as.numeric(table(factor(unlist(proteins[[sp]]), levels=a())))
  } else{
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
rownames(freqaa)[which(rownames(freqaa)=="Pauxi_pauxi")]="BRAKER2"

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
