####################################################################

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=2){
  stop("You need to provide 5 arguments: cluster, annotation, sample, nthreads, strand")
}

sample=args[1]
nthreads=as.numeric(args[2])

####################################################################

library(Rsubread)

####################################################################

path="/sps/biometr/necsulea/HelmetedCurassowGenome/"
pathRNASeq=paste(path, "data/RNASeq/Mouse/",sep="")
pathAnnot=paste(path,"data/genome_annotations/Ensembl102/Mus_musculus.GRCm38.102.chr.formatted.gtf",sep="")
pathResults=paste(path, "results/FaceBase_analysis/Mouse/",sep="")
  
####################################################################

file=paste(sample,"_sorted.bam",sep="")
pathAln=paste(pathRNASeq, file, sep="")

res.unique=featureCounts(files=pathAln, annot.ext=pathAnnot, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", countMultiMappingReads=FALSE, nthreads=nthreads)

counts=as.data.frame(res.unique$counts)
colnames(counts)=sample

annot=as.data.frame(res.unique$annotation)
rownames(annot)=as.character(annot$GeneID)

counts$Length=annot[rownames(counts), "Length"]

####################################################

write.table(counts, file=paste(pathResults, "ReadCounts_", sample, ".txt", sep=""), row.names=F, col.names=T, sep="\t")

####################################################################



