####################################################################

species=c("Basiliscus_vittatus", "Pauxi_pauxi")
reflist=list()
reflist[["Basiliscus_vittatus"]]=c("Anolis_carolinensis")
reflist[["Pauxi_pauxi"]]=c("Chicken", "Duck")

ensrelease=103
minpcid=60
assembly="MEGAHIT_RAGOUT"

col.set=c("black", "red", "royalblue", "lightseagreen", "violet",  "darkgoldenrod", "turquoise4", "yellowgreen", "indianred", "darkred","yellow", "orange", "darkorange", "darkseagreen", "darkturquoise", "deeppink",  "goldenrod4", "darkslategray", "darkslateblue", "lightblue1", "khaki2", "green", "hotpink",  "yellow3", "lightsalmon", "mediumpurple1",  "navajowhite2", "lightsteelblue2", "steelblue", "saddlebrown","peru")

pathAnnot="../../data/ensembl_annotations/"
pathFigures="../../results/figures/"

####################################################################

for(sp in species){

    pathGenome=paste("../../results/genome_assembly/",sp,"/",assembly,"/",sep="")

    seqnames=read.table(paste(pathGenome, "sequence_names.txt",sep=""), h=F, stringsAsFactors=F)
    rownames(seqnames)=seqnames[,2]

    synonyms=seqnames[,2]
    names(synonyms)=seqnames[,1]

    revsynonyms=seqnames[,1]
    names(revsynonyms)=seqnames[,2]

    chrsizes=read.table(paste(pathGenome, "GCContent_ChromosomeSize.txt",sep=""), h=T, stringsAsFactors=F)
    rownames(chrsizes)=chrsizes$Chr

    chrsizes=chrsizes[order(chrsizes$Size, decreasing=T),]

    sizes=chrsizes$Size
    names(sizes)=seqnames[chrsizes$Chr, 1]

    pathResults=paste("../../results/genome_assembly_quality/",sp,"/",assembly,"/",sep="")

    for(ref in reflist[[sp]]){

        ## gene coordinates

        genecoords=read.table(paste(pathAnnot, ref, "/GeneCoordinates_Ensembl", ensrelease,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

        genecoords$MeanPos=apply(genecoords[, c("Gene.start..bp.","Gene.end..bp.")], 1, mean)
        genecoords=genecoords[which(!duplicated(genecoords$Protein.stable.ID)),]
        rownames(genecoords)=genecoords$Protein.stable.ID

        ## tblastn results

        aln=read.table(paste(pathResults, "AlignmentStatistics_", ref, "_AllPeptides",ensrelease,"_vs_genome_sequence_minPCIdentity",minpcid,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

        aln$ProteinID=unlist(lapply(aln$ProteinID, function(x) unlist(strsplit(x, split="\\."))[1]))
        aln$GeneID=unlist(lapply(aln$GeneID, function(x) unlist(strsplit(x, split="\\."))[1]))

        aln$PropAln=aln$AlignedLength/aln$TotalLength
        aln=aln[which(aln$PropAln>0.5),]
        aln=aln[order(aln$PropAln, decreasing=T),]
        aln=aln[which(!duplicated(aln$GeneID)),] ## best aligned protein per gene

        aln$RefChr=genecoords[aln$ProteinID,"Chromosome.scaffold.name"]
        aln$RefStrand=genecoords[aln$ProteinID,"Strand"]
        aln$RefPos=genecoords[aln$ProteinID,"MeanPos"]

        ## chr colors
        nb.hits=as.numeric(table(as.factor(aln$RefChr)))
        names(nb.hits)=levels(as.factor(aln$RefChr))
        nb.hits=sort(nb.hits,decreasing=T)
        all.refchr=names(nb.hits)
        all.colors=rep("gray60", length(all.refchr))
        names(all.colors)=all.refchr

        atleast50=length(which(nb.hits>=50))
        all.colors[1:min(length(col.set), atleast50)]=col.set[1:min(length(col.set), atleast50)]

        aln=aln[which(!is.na(aln$Start)),] ## unique hits

        aln$TgPos=(aln$Start+aln$End)/2

        ## nb aligned genes per contig

        nb.genes.tg=as.numeric(table(aln$Contig))
        names(nb.genes.tg)=levels(as.factor(aln$Contig))
        nb.genes.tg=sort(nb.genes.tg, decreasing=T)

        nb.genes.ref=as.numeric(table(aln$RefChr))
        names(nb.genes.ref)=levels(as.factor(aln$RefChr))
        nb.genes.ref=sort(nb.genes.ref, decreasing=T)

        ## select colors

        sel.chr=names(nb.genes.ref)[which(nb.genes.ref>=50)]

        chr.legend=intersect(c(as.character(1:33), "Z", "W"), sel.chr)
        n=length(chr.legend)

        col.vector=col.set[1:length(chr.legend)]
        names(col.vector)=chr.legend

        nb.chr=length(chr.legend)
        mean.nb.chr=floor(nb.chr/2)

        ## select contigs with at least 50 genes with tblastn hits
        ## contigs=names(nb.genes.tg)[which(nb.genes.tg>=50)]

        ## select contigs with at least 3 Mb
        selcontigs=chrsizes$Chr[which(chrsizes$Size>=3e6)]
        contigs=revsynonyms[selcontigs]

        ## do plot

        pdf(file=paste(pathFigures, "Synteny_",sp,"_",ref,".pdf",sep=""), width=10, height=10)
        par(mar=c(3.1, 2.1,1.1, 5.1))

        minx=1
        maxx=max(aln$End)

        ypos=length(contigs):1
        smally=2/length(ypos)

        names(ypos)=contigs

        plot(1, type="n", xlim=c(minx, maxx), ylim=c(0,max(ypos)+1), axes=F, xlab="", ylab="", main="")

        for(contig in contigs){
            this.aln=aln[which(aln$Contig==contig),]

            this.size=sizes[contig]
            this.syn=synonyms[contig]

            xpos=this.aln$TgPos

            segments(minx, ypos[contig], maxx, ypos[contig], lty=2, col="gray40")

            rect(minx, ypos[contig]-smally, this.size, ypos[contig]+smally, col="white", border="gray40")

            points(this.aln$TgPos, rep(ypos[contig], nrow(this.aln)), pch=20, col=all.colors[this.aln$RefChr], cex=0.6)

            mtext(this.syn, side=4, las=2, at=ypos[contig])
        }

        prettyx=pretty(c(minx, maxx)/1e6)
        axislab=paste(prettyx, "Mb",sep="")

        axis(side=3, pos=length(contigs)+0.5, at=prettyx*1e6, labels=axislab)

        if(nb.chr>6){
            legend("bottomleft", horiz=T, legend=chr.legend[1:mean.nb.chr], pch=20, col=all.colors[1:mean.nb.chr], inset=c(0, -0.001), bty="n", cex=0.9, xpd=NA, pt.cex=2, title=paste(ref, "chromosomes"))
            legend("bottomleft", horiz=T, legend=c(chr.legend[(mean.nb.chr+1):nb.chr], "other"), pch=20, col=c(all.colors[(mean.nb.chr+1):nb.chr], "gray60"), inset=c(0, -0.05), bty="n",xpd=NA, cex=0.9, pt.cex=2)
        } else{
            legend("bottomleft", horiz=T, legend=c(chr.legend[1:nb.chr], "other"), pch=20, col=c(all.colors[1:nb.chr], "gray60"), inset=c(0, -0.001), bty="n", cex=0.9, xpd=NA, pt.cex=2, title=paste(ref, "chromosomes"))
        }

        dev.off()
    }
}

####################################################################
