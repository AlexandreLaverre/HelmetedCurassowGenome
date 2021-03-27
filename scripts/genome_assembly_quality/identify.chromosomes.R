####################################################################

pathResults="../../results/genome_assembly_quality/MEGAHIT_RAGOUT/"

####################################################################

corresp.gga=read.table(paste(pathResults, "ChromosomeCorrespondence_Chicken.txt", sep=""), h=T, stringsAsFactors=F)
corresp.apl=read.table(paste(pathResults, "ChromosomeCorrespondence_Duck.txt", sep=""), h=T, stringsAsFactors=F)

rownames(corresp.gga)=corresp.gga$Contig.Scaffold
rownames(corresp.apl)=corresp.apl$Contig.Scaffold

## select reciprocal best hits, at least 5 proteins

corresp.gga=corresp.gga[which(corresp.gga$Contig.Scaffold==corresp.gga$TargetBestHit & corresp.gga$NbSupportingProteins>=5),]
corresp.apl=corresp.apl[which(corresp.apl$Contig.Scaffold==corresp.apl$TargetBestHit & corresp.apl$NbSupportingProteins>=5),]
  
####################################################################

cs=unique(c(corresp.gga$Contig.Scaffold, corresp.apl$Contig.Scaffold))

results=data.frame("Contig.Scaffold"=cs, "Chicken"=corresp.gga[cs,"BestHit"], "Duck"=corresp.apl[cs, "BestHit"], stringsAsFactors=F)

results=results[order(results$Chicken),]

write.table(results, file=paste(pathResults, "ChromosomeCorrespondence_Chicken_Duck.txt", sep=""), row.names=F, col.names=T, sep="\t")

####################################################################

