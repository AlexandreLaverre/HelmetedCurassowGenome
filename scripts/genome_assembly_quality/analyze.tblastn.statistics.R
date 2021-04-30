#########################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

pathResults=paste(path, "results/genome_assembly_quality/",sep="")
pathAnnot=paste(path, "data/ensembl_annotations/",sep="")
pathFigures=paste(path, "results/figures/", sep="")

#########################################################################

method="MEGAHIT_RAGOUT"

ensrelease=103

#########################################################################

mat.counts.50=matrix(rep(NA,3*4), nrow=3)
rownames(mat.counts.50)=c("Chicken", "Duck", "PaintedTurtle")
colnames(mat.counts.50)=as.character(c(50, 60, 70, 80))

mat.prop.50=mat.counts.50
mat.counts.0=mat.counts.50
mat.prop.0=mat.counts.50

for(ref in c("Chicken", "Duck", "PaintedTurtle")){
  
  genecoords=read.table(paste(pathAnnot, ref, "/GeneCoordinates_Ensembl", ensrelease,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  genecoords=genecoords[which(genecoords$Protein.stable.ID!=""),]
  nbtot=length(unique(genecoords$Gene.stable.ID))
  
  
  for(minpcid in c(50, 60, 70, 80)){
    aln=read.table(paste(pathResults, method,"/AlignmentStatistics_", ref, "_AllPeptides",ensrelease,"_vs_genome_sequence_minPCIdentity",minpcid,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
    
    aln$ProteinID=unlist(lapply(aln$ProteinID, function(x) unlist(strsplit(x, split="\\."))[1]))
    aln$GeneID=unlist(lapply(aln$GeneID, function(x) unlist(strsplit(x, split="\\."))[1]))
    
    aln$PropAln=aln$AlignedLength/aln$TotalLength
    aln=aln[order(aln$PropAln, decreasing=T),]
    aln=aln[which(!duplicated(aln$GeneID)),] ## best aligned protein per gene
    
    nbaln50=length(which(aln$PropAln>0.5))
    nbaln0=length(unique(aln$GeneID))

    mat.counts.0[ref,as.character(minpcid)]=nbaln0
    mat.counts.50[ref,as.character(minpcid)]=nbaln50
    mat.prop.0[ref,as.character(minpcid)]=nbaln0/nbtot
    mat.prop.50[ref,as.character(minpcid)]=nbaln50/nbtot
  }
}

pdf(file=paste(pathFigures, "NbGenesWithTBLASTNHits_",method,".pdf",sep=""), width=8, height=5)
par(mfrow=c(1,2))
par(mar=c(4.1, 4.5, 4.1, 1.1))
barplot(mat.counts.0,beside=T, col=c("indianred", "darkred", "steelblue"), ylim=c(0, 20000), xlab="minimum % identity", ylab="nb. genes", cex.lab=1.1)
legend("topleft", c("chicken", "duck", "painted turtle"), horiz=T, fill=c("indianred", "darkred", "steelblue"), inset=c(0, -0.25), bty="n",xpd=NA, cex=1.1) 
barplot(100*mat.prop.0,beside=T, col=c("indianred", "darkred", "steelblue"), ylim=c(0, 100), xlab="minimum % identity", ylab="% genes", cex.lab=1.1)

dev.off()

#########################################################################
