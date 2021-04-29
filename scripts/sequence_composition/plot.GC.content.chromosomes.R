#########################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

pathGenomeAssembly=paste(path, "results/genome_assembly/",sep="")
pathFigures=paste(path, "results/figures/", sep="")

#########################################################################

for(method in c("MEGAHIT", "MEGAHIT_RAGOUT")){
  gc=read.table(paste(pathGenomeAssembly, method, "/GCContent_ChromosomeSize.txt", sep=""), h=T, stringsAsFactors=F)

  print(paste("average GC content", round(mean(gc$GC), digits=2)))

  pdf(file=paste(pathFigures, "GCContent_Contigs_",method,".pdf", sep=""), width=6, height=5.5)

  hist(gc$GC, xlab="contig/scaffold GC content", main=method)
  
  dev.off()
}

#########################################################################
