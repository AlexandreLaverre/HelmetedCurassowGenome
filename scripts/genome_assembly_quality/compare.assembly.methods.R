#########################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

pathGenomeAssembly=paste(path, "results/genome_assembly_tests/",sep="")
pathAnnotations=paste(path, "data/ensembl_annotations/", sep="")

#########################################################################

genomeprefix="chr4"
ensrelease="98"

#########################################################################

geneinfo=read.table(paste(pathAnnotations, "Chicken/GeneInfo_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

genes.chr4=geneinfo$stable_id[which(geneinfo$name=="4")]

kmers=c(23, 27, 33, 37, 43, 47, 53, 57, 63)

#########################################################################

coverage=list()
breakpoints=list()
goodproteins=list()
n50=list()

for(method in c("SOAPdenovo", "ABYSS")){

  coverage[[method]]=c()
  breakpoints[[method]]=c()
  goodproteins[[method]]=c()
  n50[[method]]=c()
  
  for(kmer in kmers){
    tblastn.statsfile=system(paste("ls ",pathGenomeAssembly, method, "/",genomeprefix," | grep kmer", kmer, " | grep tblastn.stats.out", sep=""), intern=T)

    if(length(tblastn.statsfile)==1 & file.exists(paste(pathGenomeAssembly, method, "/",genomeprefix, "/", tblastn.statsfile,sep=""))){
      this.tblastn.stats=read.table(paste(pathGenomeAssembly, method, "/",genomeprefix, "/", tblastn.statsfile,sep=""), h=T, stringsAsFactors=F, sep="\t")
      
      this.tblastn.stats$GeneID=unlist(lapply(this.tblastn.stats$GeneID, function(x) unlist(strsplit(x, split="\\."))[1]))
      
      this.tblastn.stats=this.tblastn.stats[which(this.tblastn.stats$GeneID%in%genes.chr4),]
      
      this.tblastn.stats$AlignedRatio=this.tblastn.stats$AlignedLength/this.tblastn.stats$TotalLength
      
      maxratio.gene=tapply(this.tblastn.stats$AlignedRatio, as.factor(this.tblastn.stats$GeneID), max)
      
      prop.good=round(length(which(maxratio.gene>=0.75))/length(maxratio.gene), digits=3)
      prop.verygood=round(length(which(maxratio.gene>=0.9))/length(maxratio.gene), digits=3)

      goodproteins[[method]]=c(goodproteins[[method]], prop.good)
    } else{
      goodproteins[[method]]=c(goodproteins[[method]], NA)
    }

    assembly.statsfile=system(paste("ls ",pathGenomeAssembly, method, "/",genomeprefix," | grep kmer", kmer, " | grep assembly.stats.out", sep=""), intern=T)

    if(length(assembly.statsfile)==1 & file.exists(paste(pathGenomeAssembly, method, "/",genomeprefix, "/", assembly.statsfile,sep=""))){
      this.ass.stats=read.table(paste(pathGenomeAssembly, method, "/",genomeprefix, "/", assembly.statsfile,sep=""), h=T, stringsAsFactors=F, sep="\t")
      
      this.n50=this.ass.stats$N50
      this.coverage=this.ass.stats$GenomeCoverage
      this.breakpoints=this.ass.stats$NbBreakpoints
      
      n50[[method]]=c(n50[[method]], this.n50)
      coverage[[method]]=c(coverage[[method]], this.coverage)
      breakpoints[[method]]=c(breakpoints[[method]], this.breakpoints)
    } else{
      n50[[method]]=c(n50[[method]], NA)
      coverage[[method]]=c(coverage[[method]], NA)
      breakpoints[[method]]=c(breakpoints[[method]], NA)
    }
  }
}

#########################################################################

## plot coverage

pdf(file="figures/AssemblyStatistics_Coverage_chr4.pdf", width=8, height=8)

m=matrix(c(coverage[["ABYSS"]], coverage[["SOAPdenovo"]]), nrow=2, byrow=T)

par(mar=c(4.1,4.1,3.1,2.1))
b=barplot(m, beside=T, ylim=c(0, 1e8), col=c("gray40", "gray90"))

xpos=apply(b,2, mean)
mtext(kmers, at=xpos, side=1, line=1)
mtext("Genome coverage", side=2, line=3)

legend("topright", c("ABYSS", "SOAPdenovo"), fill=c("gray40", "gray90"), inset=c(0.01,-0.075), xpd=NA)

dev.off()

#########################################################################

pdf(file="figures/AssemblyStatistics_N50_chr4.pdf", width=8, height=8)

m=matrix(c(n50[["ABYSS"]], n50[["SOAPdenovo"]]), nrow=2, byrow=T)

par(mar=c(4.1,4.1,3.1,2.1))
b=barplot(m, beside=T, ylim=c(0, 2300000), col=c("gray40", "gray90"))

xpos=apply(b,2, mean)
mtext(kmers, at=xpos, side=1, line=1)
mtext("N50", side=2, line=3)

legend("topleft", c("ABYSS", "SOAPdenovo"), fill=c("gray40", "gray90"), inset=c(0.01,-0.01), xpd=NA)

dev.off()

#########################################################################


pdf(file="figures/AssemblyStatistics_Breakpoints_chr4.pdf", width=8, height=8)

m=matrix(c(breakpoints[["ABYSS"]], breakpoints[["SOAPdenovo"]]), nrow=2, byrow=T)

par(mar=c(4.1,4.1,3.1,2.1))
b=barplot(m, beside=T, ylim=c(0, 120000), col=c("gray40", "gray90"))

xpos=apply(b,2, mean)
mtext(kmers, at=xpos, side=1, line=1)
mtext("nb. breakpoints", side=2, line=3)

legend("topright", c("ABYSS", "SOAPdenovo"), fill=c("gray40", "gray90"), inset=c(0.01,0.01), xpd=NA)

dev.off()

#########################################################################

pdf(file="figures/AssemblyStatistics_ReconstructedProteins_chr4.pdf", width=8, height=8)

m=100*matrix(c(goodproteins[["ABYSS"]], goodproteins[["SOAPdenovo"]]), nrow=2, byrow=T)

par(mar=c(4.1,4.1,3.1,2.1))
b=barplot(m, beside=T, ylim=c(0, 100), col=c("gray40", "gray90"))

xpos=apply(b,2, mean)
mtext(kmers, at=xpos, side=1, line=1)
mtext("% reconstructed proteins", side=2, line=3)

legend("topright", c("ABYSS", "SOAPdenovo"), fill=c("gray40", "gray90"), inset=c(0.01,-0.075), xpd=NA)

dev.off()

#########################################################################
