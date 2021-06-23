#######################################################################

pathGeMoMa="../../results/genome_annotation/MEGAHIT_RAGOUT/GeMoMa/combined/"
pathOrthoFinder=paste(pathGeMoMa, "OrthoFinder_May18/Results_May18/", sep="")
pathFigures="../../results/figures/"

#######################################################################

library(seqinr)
library(ade4)

#######################################################################

proteins=read.fasta(paste(pathGeMoMa, "primary_transcripts/filtered_predictions_formatted.faa", sep=""), seqtype="AA")

protlen=unlist(lapply(proteins, length))

#######################################################################

orthogroups=read.table(paste(pathOrthoFinder, "Phylogenetic_Hierarchical_Orthogroups/N0.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

ortho=setdiff(orthogroups$Pauxi_pauxi, "")
ortho=unlist(lapply(ortho, function(x) unlist(strsplit(x, split=", "))))

notortho=setdiff(names(proteins), ortho)

#######################################################################

## protein length

pdf(paste(pathFigures, "ProteinLength_GeMoMa_Orthogroups.pdf",sep=""))

do=density(protlen[ortho], bw=38)
dn=density(protlen[notortho], bw=38)
xlim=range(c(do$x, dn$x))
ylim=range(c(do$y, dn$y))

xlim[2]=2000

par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(do$x, do$y, xlim=xlim, ylim=ylim, type="l", xlab="protein length", ylab="density", col="black")
lines(dn$x, dn$y, col="red")

legend("topright", inset=0.01, col=c("black", "red"), legend=c("in orthogroups", "not in orthogroups"), lty=1)

dev.off()

#######################################################################

## select proteins with at least 100 aa
long=names(proteins)[which(protlen>=100)]

#######################################################################

ovrep=read.table(paste(pathGeMoMa, "overlap_repeats.txt", sep=""), h=T, stringsAsFactors=F)
ovrep$PercentageRepeats=100*ovrep$RepeatOverlapLength/ovrep$TranscriptLength

ovrepgene=tapply(ovrep$PercentageRepeats, ovrep$GeneID, max)

#######################################################################

## protein length

pdf(paste(pathFigures, "OverlapRepeats_GeMoMa_Orthogroups.pdf",sep=""))

do=density(ovrepgene[ortho], bw=0.5, na.rm=T)
dn=density(ovrepgene[notortho], bw=0.5, na.rm=T)
xlim=range(c(do$x, dn$x))
ylim=range(c(do$y, dn$y))

par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(do$x, do$y, xlim=xlim, ylim=ylim, type="l", xlab="% repeats", ylab="density", col="black")
lines(dn$x, dn$y, col="red")

legend("topright", inset=0.01, col=c("black", "red"), legend=c("in orthogroups", "not in orthogroups"), lty=1)

dev.off()

#######################################################################

## read protein sequences for all species

files=system(paste("ls ",pathGeMoMa,"OrthoFinder_May18/ | grep fa$", sep=""), intern=T)
species=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))
species=setdiff(species, "Pauxi_pauxi")

################################################################################

freqaa=list()

for(sp in species){
  print(sp)
  
  this.proteins=read.fasta(paste(pathGeMoMa,"OrthoFinder_May18/", sp, ".fa", sep=""), seqtype="AA")

  freqaa[[sp]]=as.numeric(table(factor(unlist(this.proteins), levels=a())))
}

################################################################################

## add hocco

freqaa[["all"]]=as.numeric(table(factor(unlist(proteins), levels=a())))
freqaa[["ortho"]]=as.numeric(table(factor(unlist(proteins[ortho]), levels=a())))
freqaa[["long"]]=as.numeric(table(factor(unlist(proteins[long]), levels=a())))

################################################################################

table.freqaa=t(as.data.frame(freqaa))
colnames(table.freqaa)=a()
rownames(table.freqaa)=names(freqaa)

################################################################################

## COA, all proteins

afc <- dudi.coa(table.freqaa, scann = FALSE, nf = 5)

pdf(paste(pathFigures, "COA_AminoAcidComposition_GeMoMa_OrthoGroups.pdf",sep=""), width=6, height=6)

par(mar=c(4.5, 5.1, 2.1, 1.1))
plot(afc$li[,1],afc$li[,2], pch = 19, main = "correspondence analysis, first factorial map", xlab = "", ylab ="",las = 1)

mtext(paste0("F1 (", round(100*afc$eig[1]/sum(afc$eig)), "% explained variance)"), side=1, line=3)
mtext(paste0("F2 (", round(100*afc$eig[2]/sum(afc$eig)), "% explained variance)"), side=2, line=3.5)

smallx=diff(range(afc$li[,1]))/100
smally=diff(range(afc$li[,2]))/100

points(afc$li["all",1],afc$li["all",2], pch=19, col="red")
points(afc$li["ortho",1],afc$li["ortho",2], pch=19, col="steelblue")
points(afc$li["long",1],afc$li["long",2], pch=19, col="seagreen")


legend("bottomleft", legend=c("GeMoMa, all proteins", "GeMoMa, in orthogroups", "GeMoMa, min 100aa"), pch=20, col=c("red", "steelblue", "seagreen"), inset=0.01)


dev.off()

################################################################################
