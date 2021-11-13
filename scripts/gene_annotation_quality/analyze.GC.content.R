#############################################################################

pathAnnot="../../results/genome_annotation/"
pathFigures="../../results/figures/"

assembly="MEGAHIT_RAGOUT"
gemoma="GeMoMa"

##############################################################################

library(seqinr)
library(ade4)

#############################################################################

cds=read.fasta(paste(pathAnnot, assembly, "/GeMoMa/combined/filtered_predictions.cds.fa", sep=""), seqtype="DNA")

cds.len=unlist(lapply(cds, length))

cds=cds[which(cds.len>100)]

#############################################################################

cds.gc=unlist(lapply(cds, GC))

#############################################################################

pdf(file=paste(pathFigures, "GC_content_GeMoMa.pdf",sep=""))

hist(100*cds.gc, breaks=50, main="", xlab="GC content")
dev.off()

#############################################################################

