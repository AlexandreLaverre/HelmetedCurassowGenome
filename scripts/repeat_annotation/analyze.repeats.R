#########################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

pathRepeats=paste(path, "results/repeats/",sep="")
pathFigures=paste(path, "results/figures/", sep="")

#########################################################################

method="MEGAHIT_RAGOUT"

dr=read.table(paste(pathRepeats, method, "/RepeatMasker/Comparison_Dfam_RepeatModeler.txt", sep=""), h=T, stringsAsFactors=F)
rd=read.table(paste(pathRepeats, method,  "/RepeatMasker/Comparison_RepeatModeler_Dfam.txt", sep=""), h=T, stringsAsFactors=F) 
combined=read.table(paste(pathRepeats,  method, "/RepeatMasker/CombinedRepeatAnnotations.txt", sep=""), h=T, stringsAsFactors=F)

dr$Family=unlist(lapply(dr$Type, function(x) unlist(strsplit(x, split="\\/"))[1]))
rd$Family=unlist(lapply(rd$Type, function(x) unlist(strsplit(x, split="\\/"))[1]))

dr$Family[which(dr$Family=="DNA?")]="DNA"

#########################################################################

all.types=unique(c(dr$Type, rd$Type))
all.families=unique(c(dr$Family, rd$Family))

top=c("DNA", "LINE", "SINE", "LTR")
simple=c("Simple_repeat", "Low_complexity")
RNA=c("scRNA", "snRNA", "sprRNA", "rRNA", "tRNA")
bottom=c("Unknown")
  
all.families=c(top, simple, RNA, setdiff(all.families, c(top, simple, RNA,bottom)), bottom)

#########################################################################

library(RColorBrewer)

n=length(all.families)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col.vector=sample(col.vector, n)
names(col.vector)=all.families

####################################################################

nb.dfam=table(factor(dr$Family, levels=all.families))
nb.rm=table(factor(rd$Family, levels=all.families))

m=matrix(c(nb.dfam, nb.rm), nrow=2, byrow=T)

pdf(file=paste(pathFigures, "RepeatClasses_",method,".pdf",sep=""), width=4.5, height=6)
par(mar=c(3.1, 2.1,1.1, 8.1))

b=barplot(t(m), beside=F, col=col.vector)
mtext(c("Dfam", "RepeatModeler"), side=1, at=b, line=0.5)

legend("topright", legend=all.families[n:1], fill=col.vector[n:1], bty="n", inset=c(-0.65, 0), xpd=NA)

dev.off()

####################################################################

dr$Size=dr$End-dr$Start+1
rd$Size=rd$End-rd$Start+1

combined$Size=combined$End-combined$Start+1

print(paste("Dfam", sum(dr$Size/1e6)))
print(paste("RepeatModeler", sum(rd$Size/1e6)))
print(paste("combined", sum(combined$Size/1e6)))

####################################################################
