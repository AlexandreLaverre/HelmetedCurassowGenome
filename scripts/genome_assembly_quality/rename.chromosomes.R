##############################################################

pathResults="../../results/genome_assembly/MEGAHIT_RAGOUT/"

##############################################################

scaffolds=system(paste("grep \">\" ",pathResults, "Pauxi_pauxi_scaffolds.fasta", sep=""), intern=T)
contigs=system(paste("grep \">\" ",pathResults, "Pauxi_pauxi_unplaced.fasta", sep=""), intern=T)

scaffolds=unlist(lapply(scaffolds, function(x) substr(x, 2, nchar(x))))
contigs=unlist(lapply(contigs, function(x) substr(x, 2, nchar(x))))

##############################################################

s=readLines(paste(pathResults, "assembly.stats.out", sep=""))

id=unlist(lapply(s, function(x) {y=unlist(strsplit(x, split=":")); return(paste(y[1:(length(y)-1)], collapse=":"))}))
size=as.numeric(unlist(lapply(s, function(x) {y=unlist(strsplit(x, split=":")); return(y[length(y)])})))

stats=data.frame("ID"=id, "Size"=size, stringsAsFactors=F)
stats=stats[4:dim(stats)[1],]

stats$Type=rep(NA, dim(stats)[1])
stats$Type[which(stats$ID%in%scaffolds)]="Scaffold"
stats$Type[which(stats$ID%in%contigs)]="Contig"

stats=stats[order(stats$Size, decreasing=T),]
stats$NewName=rep(NA, dim(stats)[1])

stats$NewName[which(stats$Type=="Scaffold")]=paste("Scaffold", 1:length(which(stats$Type=="Scaffold")), sep="")
stats$NewName[which(stats$Type=="Contig")]=paste("Contig", 1:length(which(stats$Type=="Contig")), sep="")

##############################################################

## output

write.table(stats[,c("ID", "NewName")], file=paste(pathResults, "sequence_names.txt",sep=""), sep="\t", row.names=F, col.names=F, quote=F)

##############################################################
