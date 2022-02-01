#####################################################################

pathResults="../../results/coding_gene_evolution/"
pathPhenotypes="../../data/MGI_alleles/"
pathExp="../../results/FaceBase_analysis/Mouse/"
pathAnnot="../../data/genome_annotations/Ensembl103/"
pathFigures="../../results/figures/"

#####################################################################

mgi=read.table(paste(pathPhenotypes, "MGIalleleQuery_20220131_042202.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"", comment.char="") 

mgi$Gene=toupper(unlist(lapply(mgi$Allele.Symbol, function(x) unlist(strsplit(x, split="<"))[1])))

#####################################################################

files=system(paste("ls ",pathResults, "data_for_pelican/",sep=""), intern=T)
genes=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))
genes=unlist(lapply(genes, function(x) unlist(strsplit(x, split="_"))[3]))
genes=setdiff(genes, NA)

#####################################################################

tested.mgi.genes=intersect(mgi$Gene, genes)

#####################################################################

best_sites=read.table(paste(pathResults, "pelican_output_by_category/best_sites.tsv",sep=""), h=T, stringsAsFactors=F)
best_sites$Gene=unlist(lapply(best_sites$alignment,function(x) unlist(strsplit(x, split="_"))[3]))

best_sites=best_sites[which(best_sites$Gene%in%genes),]

#####################################################################

best_sites=best_sites[which(best_sites$aa_model_pval<0.05),]

signif.genes=unique(best_sites$Gene)

## writeLines(unique(signif.genes), con="signif_genes_0.05.txt")

nbsites=as.numeric(table(best_sites$Gene))
names(nbsites)=levels(as.factor(best_sites$Gene))

#####################################################################

gene_score=tapply(best_sites$aa_model_pval, as.factor(best_sites$Gene), function(x) sum(-log10(x)))

#####################################################################

tpm=read.table(paste(pathExp,  "TPM.txt", sep=""), h=T,stringsAsFactors=F)

meantpm=apply(tpm,1,mean)
meantpm.fnp=apply(tpm[,grep("FNP", colnames(tpm))],1,mean)

#####################################################################

names=read.table(paste(pathAnnot, "Mus_musculus.genenames.txt",sep=""),h=F, stringsAsFactors=F)
rownames(names)=names[,1]
names=names[which(names[,1]%in%rownames(tpm)),]

names[,2]=toupper(names[,2])
names=names[which(names[,2]%in%genes),]

dupli=names$V2[which(duplicated(names$V2))]
names=names[which(!names$V2%in%dupli),]

rownames(names)=names$V2

#####################################################################

exp.all=tpm[names$V1,]
rownames(exp.all)=names$V2

exp.signif=tpm[names$V1[which(names$V2%in%best_sites$Gene)],]
rownames(exp.signif)=names$V2[which(names$V2%in%best_sites$Gene)]

#####################################################################
  
tiss=unlist(lapply(colnames(exp.all), function(x) paste(unlist(strsplit(x, split="_"))[1:2], collapse="_")))

maxtiss.signif=apply(exp.signif, 1, function(x) tiss[which.max(x)])
maxexp.signif=apply(exp.signif, 1, max)
maxtiss.signif[which(maxexp.signif<1)]=NA
names(maxtiss.signif)=rownames(exp.signif)

maxtiss.all=apply(exp.all, 1, function(x) tiss[which.max(x)])
maxexp.all=apply(exp.all, 1, max)
maxtiss.all[which(maxexp.all<1)]=NA
names(maxtiss.all)=rownames(exp.all)

#####################################################################

best_sites$maxexp=maxexp.all[best_sites$Gene]
best_sites$maxtiss=maxtiss.all[best_sites$Gene]

#####################################################################

exp.samples = unlist(lapply(colnames(exp.all), function(x) paste(unlist(strsplit(x, split="_"))[-3], collapse="_")))

meanexp.all=t(apply(exp.all, 1, function(x) tapply(as.numeric(x), as.factor(exp.samples), mean)))

#####################################################################

compute.tau <- function(exp){
  if(max(exp)==0){
    return(NA)
  }
  
  n=length(exp)
  newexp=exp/max(exp)

  tau=sum(1-newexp)/(n-1)

  return(tau)
}

tau=apply(meanexp.all,1, compute.tau)

#####################################################################

best_sites$tau=tau[best_sites$Gene]

maxsample.all=apply(meanexp.all, 1, function(x) colnames(meanexp.all)[which.max(x)])

maxexp.avg.all=apply(meanexp.all,1, mean)

best_sites$maxsample=maxsample.all[best_sites$Gene]

#####################################################################

ciliopathy=toupper(readLines(paste(pathResults, "ciliophaty_genes.txt",sep="")))

corematrisome=toupper(readLines(paste(pathResults,"core_matrisome.txt",sep="")))
assocmatrisome=toupper(readLines(paste(pathResults,"matrisome_associated.txt",sep="")))
  
#####################################################################

results=data.frame("GeneName"=signif.genes, "NbSignifSites"=nbsites[signif.genes], "GeneScore"=gene_score[signif.genes], stringsAsFactors=F)


other.genes=setdiff(genes, signif.genes)

results.other=data.frame("GeneName"=other.genes, "NbSignifSites"=rep(0, length(other.genes)), "GeneScore"=rep(0, length(other.genes)), stringsAsFactors=F)

results=rbind(results, results.other, stringsAsFactors=F)

results=results[order(results$GeneScore, decreasing=T),]

#####################################################################

results$MGICraniofacial=results$GeneName%in%tested.mgi.genes
results$Ciliopathy=results$GeneName%in%ciliopathy

results$CoreMatrisome=results$GeneName%in%corematrisome
results$MatrisomeAssociated=results$GeneName%in%assocmatrisome

results$MaxTPM_MouseDevelopingFace=maxexp.avg.all[results$GeneName]
results$MaxSample_MouseDevelopingFace=maxsample.all[results$GeneName]

#####################################################################

write.table(results, file=paste(pathResults, "pelican_output_by_category/results_with_gene_info.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

#####################################################################

pdf(file=paste(pathFigures, "Pelican_GeneScore.pdf",sep=""), width=6, height=4)

par(mar=c(4.1, 4.5, 2.1, 1.1))
barplot(results$GeneScore[which(results$NbSignifSites>0)], ylab="gene score")
dev.off()

#####################################################################

pdf(file=paste(pathFigures, "Pelican_GeneInfo_2classes.pdf",sep=""), width=7, height=5)
par(mar=c(4.1, 4.5, 2.1, 1.1))

p1=length(which(results$MGICraniofacial & results$NbSignifSites>0))/length(which(results$NbSignifSites>0))
p2=length(which(results$MGICraniofacial & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))

p3=length(which(results$Ciliopathy & results$NbSignifSites>0))/length(which(results$NbSignifSites>0))
p4=length(which(results$Ciliopathy & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))

p5=length(which(results$CoreMatrisome & results$NbSignifSites>0))/length(which(results$NbSignifSites>0))
p6=length(which(results$CoreMatrisome & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))


p7=length(which(results$MatrisomeAssociated & results$NbSignifSites>0))/length(which(results$NbSignifSites>0))
p8=length(which(results$MatrisomeAssociated & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))


b=barplot(100*c(p1,p2, p3, p4, p5, p6, p7, p8), space=c(0.1, rep(c(0.25, 2),3), 0.25), col=c("darkred", "gray40"), ylab="% genes")

mtext("MGI\ncraniofacial", side=1, line=1, at=mean(b[1:2]))
mtext("ciliopathy", side=1, line=1, at=mean(b[3:4]))
mtext("core matrisome", side=1, line=1, at=mean(b[5:6]))
mtext("matrisome-associated", side=1, line=1, at=mean(b[7:8]))

legend("topright", fill=c("darkred", "gray40"), legend=c("significant", "not significant"), inset=0.01)

dev.off()

#####################################################################

pdf(file=paste(pathFigures, "Pelican_GeneInfo_3classes.pdf",sep=""), width=7, height=5)

par(mar=c(4.1, 4.5, 2.1, 1.1))

p1=length(which(results$MGICraniofacial & results$GeneScore>2))/length(which(results$GeneScore>2))
p2=length(which(results$MGICraniofacial & results$NbSignifSites>0 & results$GeneScore<=2))/length(which(results$NbSignifSites>0 & results$GeneScore<=2))
p3=length(which(results$MGICraniofacial & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))


p4=length(which(results$Ciliopathy & results$GeneScore>2))/length(which(results$GeneScore>2))
p5=length(which(results$Ciliopathy & results$NbSignifSites>0 & results$GeneScore<=2))/length(which(results$NbSignifSites>0 & results$GeneScore<=2))
p6=length(which(results$Ciliopathy & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))


p7=length(which(results$CoreMatrisome & results$GeneScore>2))/length(which(results$GeneScore>2))
p8=length(which(results$CoreMatrisome & results$NbSignifSites>0 & results$GeneScore<=2))/length(which(results$NbSignifSites>0 & results$GeneScore<=2))
p9=length(which(results$CoreMatrisome & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))


p10=length(which(results$MatrisomeAssociated & results$GeneScore>2))/length(which(results$GeneScore>2))
p11=length(which(results$MatrisomeAssociated & results$NbSignifSites>0 & results$GeneScore<=2))/length(which(results$NbSignifSites>0 & results$GeneScore<=2))
p12=length(which(results$MatrisomeAssociated & results$NbSignifSites==0))/length(which(results$NbSignifSites==0))



b=barplot(100*c(p1,p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12), space=c(0.1, rep(c(0.25, 0.25, 2),3), 0.25, 0.25), col=c("darkred", "brown", "gray40"), ylab="% genes")

mtext("MGI\ncraniofacial", side=1, line=1, at=mean(b[1:3]))
mtext("ciliopathy", side=1, line=1, at=mean(b[4:6]))
mtext("core matrisome", side=1, line=1, at=mean(b[7:9]))
mtext("matrisome-associated", side=1, line=1, at=mean(b[10:12]))

legend("topright", fill=c("darkred", "brown", "gray40"), legend=c("significant, score > 2", "significant, score < 2",  "not significant"), inset=0.01)

dev.off()

#####################################################################
