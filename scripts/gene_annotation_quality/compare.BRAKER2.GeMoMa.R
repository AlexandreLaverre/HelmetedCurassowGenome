################################################################################

pathAnnot="../../results/genome_annotation/"
pathFigures="../../results/figures/"

assembly="MEGAHIT_RAGOUT"
gemoma="GeMoMa"

################################################################################

library(seqinr)
library(ade4)

################################################################################

## read protein sequences for all species, gemoma

files=system(paste("ls ",pathAnnot, assembly, "/GeMoMa/combined/OrthoFinder_filtered/ | grep fa$ ", sep=""), intern=T)
species=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))

################################################################################

proteins=list()
freqaa=list()

for(sp in species){
  print(sp)
  
  if(sp=="Pauxi_pauxi"){
    proteins[["all"]]=read.fasta(paste(pathAnnot, assembly, "/GeMoMa/combined/filtered_predictions.faa", sep=""), seqtype="AA")
    freqaa[["all"]]=as.numeric(table(factor(unlist(proteins[["all"]]), levels=a())))

    proteins[["first_filter"]]=read.fasta(paste(pathAnnot, assembly, "/GeMoMa/combined/filtered_predictions_orthogroups_minLength100_maxFractionRepeats0.5.faa", sep=""), seqtype="AA")
    freqaa[["first_filter"]]=as.numeric(table(factor(unlist(proteins[["first_filter"]]), levels=a())))
    
  } else{
    proteins[[sp]]=read.fasta(paste(pathAnnot, assembly, "/GeMoMa/combined/OrthoFinder_filtered/", sp, ".fa", sep=""), seqtype="AA")
    freqaa[[sp]]=as.numeric(table(factor(unlist(proteins[[sp]]), levels=a())))
  }  
}

################################################################################

## add BRAKER2
print("BRAKER2")

proteins[["BRAKER2"]]=read.fasta(paste(pathAnnot, assembly, "/BRAKER_Ensembl103/braker.faa", sep=""), seqtype="AA")
freqaa[["BRAKER2"]]=as.numeric(table(factor(unlist(proteins[["BRAKER2"]]), levels=a())))

################################################################################

freqaa=t(as.data.frame(freqaa))
colnames(freqaa)=a()
rownames(freqaa)=names(proteins)

################################################################################
                 
print("AFC")

afc <- dudi.coa(freqaa, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_GeMoMa_BRAKER2.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc$li[,1],afc$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc$eig[1]/sum(afc$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc$eig[2]/sum(afc$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc$li[,1]))/100
smally=diff(range(afc$li[,2]))/100


points(afc$li["BRAKER2",1],afc$li["BRAKER2",2], pch=20, col="seagreen")
points(afc$li["all",1],afc$li["all",2], pch=20, col="indianred")
points(afc$li["first_filter",1],afc$li["first_filter",2], pch=20, col="blue")

legend("bottomleft", legend=c("GeMoMa all", "GeMoMa filtered", "BRAKER2"), pch=20, col=c("indianred", "blue", "seagreen"), inset=0.01)

dev.off()

################################################################################



################################################################################
