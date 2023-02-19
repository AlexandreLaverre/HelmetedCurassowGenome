################################################################

library(seqinr)

################################################################

pathAnnot="../../results/genome_annotation/"
pathResults="../../results/genome_assembly/"

################################################################

species=c("Pauxi_pauxi", "Basiliscus_vittatus")

################################################################

for(i in 1:length(species)){
    this.sp=species[i]

    ## gene coords
    gtf.sp=read.table(paste(pathAnnot, this.sp, "/MEGAHIT_RAGOUT/GeMoMa/combined/filtered_GeMoMa_annotations.gtf",sep=""), quote="", sep="\t", h=F)
    gtf.sp=gtf.sp[which(gtf.sp[,3]=="transcript"),]
    info.sp=lapply(gtf.sp[,9], function(x) unlist(strsplit(x, split=";")))
    geneid.sp=unlist(lapply(info.sp, function(x) grep("gene_id", x, value=T)))
    geneid.sp=unlist(lapply(geneid.sp, function(x) unlist(strsplit(x,split="\""))[2]))

    gtf.sp$GeneID=geneid.sp
    chr.gene=tapply(gtf.sp[,1], as.factor(gtf.sp$GeneID), function(x) unique(x)[1])
    start.gene=tapply(gtf.sp[,4], as.factor(gtf.sp$GeneID), min)
    end.gene=tapply(gtf.sp[,5], as.factor(gtf.sp$GeneID), max)
    strand.gene=tapply(gtf.sp[,7], as.factor(gtf.sp$GeneID), function(x) unique(x)[1])

    coords.sp=data.frame("GeneID"=levels(as.factor(gtf.sp$GeneID)), "Chr"=chr.gene, "Start"=start.gene, "End"=end.gene, "Strand"=strand.gene)
    rownames(coords.sp)=coords.sp$GeneID

    ## CDS
    cds=read.fasta(paste(pathAnnot, this.sp, "/MEGAHIT_RAGOUT/GeMoMa/combined/filtered_GeMoMa_annotations.cds.fa", sep=""), seqtype="DNA", forceDNAtolower=TRUE)
    genes=unlist(lapply(names(cds), function(x) paste(unlist(strsplit(x, split="_"))[1:2], collapse="_")))

    ## GC3
    cds.gc3=unlist(lapply(cds, function(x) GC3(x, forceToLower=FALSE)))
    gene.gc3=tapply(cds.gc3, as.factor(genes), mean)
    names(gene.gc3)=levels(as.factor(genes))

    df.gc3=coords.sp
    df.gc3$GC3=gene.gc3[coords.sp$GeneID]

    mean.gc3=tapply(df.gc3$GC3, as.factor(df.gc3$Chr), mean)

    df.gc3.chromo=data.frame("Chr"=levels(as.factor(df.gc3$Chr)), "MeanGC3"=mean.gc3)

    write.table(df.gc3.chromo, file=paste(pathResults, this.sp, "/MEGAHIT_RAGOUT/MeanGC3Content.txt",sep=""), row.names=F, col.names=T, quote=F)
}

################################################################
