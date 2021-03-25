#########################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

pathGenomeAssembly=paste(path, "results/genome_assembly/",sep="")
pathFigures=paste(path, "results/figures/", sep="")

#########################################################################

global.stats=list()
chr.stats=list()

for(method in c("MEGAHIT", "MEGAHIT_RAGOUT")){
  stats=read.table(paste(pathGenomeAssembly, method, "/assembly.stats.out", sep=""), h=F, stringsAsFactors=F)$V1
  id=unlist(lapply(stats, function(x) {y=unlist(strsplit(x, split=":")); return(paste(y[1:(length(y)-1)], collapse=":"))}))
  size=as.numeric(unlist(lapply(stats, function(x) {y=unlist(strsplit(x, split=":")); return(y[length(y)])})))

  stats=data.frame("ID"=id, "Size"=size, stringsAsFactors=F)

  global.stats[[method]]=stats[1:3,]
  
  this.chr.stats=stats[4:dim(stats)[1],]
  this.chr.stats=this.chr.stats[order(this.chr.stats$Size, decreasing=T),]

  chr.stats[[method]]=this.chr.stats
}

#########################################################################

pdf(file=paste(pathFigures, "AssemblyStatistics_MEGAHIT_vs_MEGAHIT_RAGOUT.pdf", sep=""), width=10, height=4.5)

par(mfrow=c(1,2))

par(mar=c(2.1, 4.1, 2.1, 1.1))

this.chr.stats=chr.stats[["MEGAHIT"]]
barplot(this.chr.stats$Size[1:50]/1e3, xlab="", ylab="contig/scaffold size (kb)",  main="MEGAHIT")
mtext("50 largest contigs/scaffolds", line=1, side=1)

this.chr.stats=chr.stats[["MEGAHIT_RAGOUT"]]
barplot(this.chr.stats$Size[1:50]/1e3, xlab="", ylab="contig/scaffold size (kb)",  main="MEGAHIT_RAGOUT")
mtext("50 largest contigs/scaffolds", line=1, side=1)

dev.off()

#########################################################################
