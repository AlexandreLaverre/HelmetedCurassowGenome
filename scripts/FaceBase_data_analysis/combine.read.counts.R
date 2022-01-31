####################################################################

path="/sps/biometr/necsulea/HelmetedCurassowGenome/"
pathRNASeq=paste(path, "data/RNASeq/Mouse/",sep="")
pathResults=paste(path, "results/FaceBase_analysis/Mouse/",sep="")
  
####################################################################

files=system(paste("ls ",pathResults," | grep ReadCounts",sep=""), intern=T)

samples=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\.txt"))[1]))
samples=unlist(lapply(samples, function(x) {y=unlist(strsplit(x, split="_")); paste(y[-1],collapse="_")}))

####################################################################

genes=c()
exp=list()
length=c()

for(sample in samples){
  rc=read.table(paste(pathResults, "ReadCounts_",sample,".txt",sep=""), h=T, stringsAsFactors=F)

  if(length(genes)==0){
    genes=rownames(rc)
    length=rc$Length
    rc=rc[genes,1]
  } else{
    rc=rc[genes,1]
  }

  exp[[sample]]=rc
}

exp=as.data.frame(exp)
M=apply(exp[,samples],2,sum)/1e6
k=length/1e3
rpk=exp/k
rpkm=t(t(exp/k)/M)
tpm=apply(rpk,2, function(x) x/(sum(x)/1e6))

rownames(rpkm)=genes
rownames(tpm)=genes


write.table(rpkm, paste(pathResults, "RPKM.txt", sep=""), row.names=T, col.names=T, sep="\t")
write.table(tpm, paste(pathResults, "TPM.txt", sep=""), row.names=T, col.names=T, sep="\t")

####################################################################
