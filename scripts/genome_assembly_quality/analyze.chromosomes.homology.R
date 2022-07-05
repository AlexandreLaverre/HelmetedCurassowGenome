####################################################################

species=c("Basiliscus_vittatus", "Pauxi_pauxi")
reflist=list()
reflist[["Basiliscus_vittatus"]]=c("Anolis_carolinensis")
reflist[["Pauxi_pauxi"]]=c("Chicken", "Duck")

ensrelease=103
minpcid=60
assembly="MEGAHIT_RAGOUT"

col.set=c("black", "red", "lightblue1", "lightseagreen", "violet",  "gray60", "gray30", "indianred", "darkred", "darkgoldenrod", "yellow", "orange", "darkorange", "darkseagreen", "darkturquoise", "deeppink",  "goldenrod4", "darkslategray", "darkslateblue", "khaki2", "green", "hotpink",  "yellow3", "lightsalmon", "mediumpurple1",  "navajowhite2", "lightsteelblue2", "steelblue", "saddlebrown","peru", "turquoise4", "yellowgreen")

pathAnnot="../../data/ensembl_annotations/"
pathFigures="../../results/figures/"

####################################################################

for(sp in species){
    for(ref in reflist[[sp]]){

        pathResults=paste("../../results/genome_assembly_quality/",sp,"/",assembly,"/",sep="")


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

        contigs=names(nb.genes.tg)[which(nb.genes.tg>=50)]

        ## do plot

        pdf(file=paste(pathFigures, "Synteny_",sp,"_",ref,".pdf",sep=""), width=10, height=10)
        par(mar=c(3.1, 2.1,1.1, 2.1))

        minx=1
        maxx=max(aln$End)

        ypos=length(contigs):1
        names(ypos)=contigs

        plot(1, type="n", xlim=c(minx, maxx), ylim=c(0,max(ypos)+1), axes=F, xlab="", ylab="", main="")

        for(contig in contigs){
            this.aln=aln[which(aln$Contig==contig),]

            xpos=this.aln$TgPos

            segments(minx, ypos[contig], maxx, ypos[contig], lty=2, col="gray40")

            points(this.aln$TgPos, rep(ypos[contig], nrow(this.aln)), pch=20, col=col.vector[this.aln$RefChr], cex=1.1)
        }

        legend("bottomleft", horiz=T, legend=chr.legend[1:mean.nb.chr], pch=20, col=col.vector[1:mean.nb.chr], inset=c(0, -0.001), bty="n", cex=0.9, xpd=NA, pt.cex=2)
        legend("bottomleft", horiz=T, legend=chr.legend[(mean.nb.chr+1):nb.chr], pch=20, col=col.vector[(mean.nb.chr+1):nb.chr], inset=c(0, -0.05), bty="n",xpd=NA, cex=0.9, pt.cex=2)


        dev.off()
    }
}

####################################################################
