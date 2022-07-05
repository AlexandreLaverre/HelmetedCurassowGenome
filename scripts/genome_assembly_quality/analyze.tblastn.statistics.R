#########################################################################

species=c("Basiliscus_vittatus", "Pauxi_pauxi")

reflist=list()
reflist[["Basiliscus_vittatus"]]=c("Anolis_carolinensis")
reflist[["Pauxi_pauxi"]]=c("Chicken", "Duck")

ensrelease=103
minpcid=50

prefixes=c("genome_sequence", "final.contigs")
names(prefixes)=c("MEGAHIT_RAGOUT", "MEGAHIT")

pathAnnot="../../data/ensembl_annotations/"
pathFigures="../../results/figures/"
pathResults="../../results/genome_assembly_quality/"

####################################################################

for(sp in species){

    for(ref in reflist[[sp]]){
        genecoords=read.table(paste(pathAnnot, ref, "/GeneCoordinates_Ensembl", ensrelease,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

        genecoords=genecoords[which(genecoords$Protein.stable.ID!=""),]
        nbtot=length(unique(genecoords$Gene.stable.ID))

        prop.aln.list=list()

        for(assembly in c("MEGAHIT_RAGOUT", "MEGAHIT")){

            prop.aln.list[[assembly]]=list()

            aln=read.table(paste(pathResults, sp, "/", assembly,"/AlignmentStatistics_", ref, "_AllPeptides",ensrelease,"_vs_",prefixes[assembly],"_minPCIdentity",minpcid,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

            aln$ProteinID=unlist(lapply(aln$ProteinID, function(x) unlist(strsplit(x, split="\\."))[1]))
            aln$GeneID=unlist(lapply(aln$GeneID, function(x) unlist(strsplit(x, split="\\."))[1]))

            aln$PropAln=aln$AlignedLength/aln$TotalLength
            aln=aln[order(aln$PropAln, decreasing=T),]
            aln=aln[which(!duplicated(aln$GeneID)),] ## best aligned protein per gene

            prop.aln.list[[assembly]]=aln$PropAln

            prop.wellaligned=round(100*length(which(aln$PropAln>0.5))/nbtot, digits=1)

            print(paste(prop.wellaligned,"% proteins are well aligned for ",sp, " vs ",ref," ",assembly,sep=""))
        }

        ## plot density of proportion of aligned sequence


        prop.aln.megahit=100*prop.aln.list[["MEGAHIT"]]
        prop.aln.ragout=100*prop.aln.list[["MEGAHIT_RAGOUT"]]

        dm=density(prop.aln.megahit, bw=1.5)
        dr=density(prop.aln.ragout, bw=1.5)

        xlim=range(c(dm$x, dr$x))
        ylim=range(c(dm$y, dr$y))

        pdf(file=paste(pathFigures, "TBLASTNHits_PCAlignedSequence_MinPCIdentity",minpcid,"_",sp,"_",ref,".pdf",sep=""), width=6, height=6)
        plot(dm$x, dm$y, col="black", xlim=xlim, ylim=ylim, xlab="% aligned protein sequence", ylab="density", type="l")
        lines(dr$x, dr$y, col="red")
        legend("topleft", lty=1, col=c("black", "red"), legend=c("MEGAHIT", "MEGAHIT_RAGOUT"), inset=0.01, cex=0.9)

        mtext(paste(sp, " vs ", ref, ", min pc identity ", minpcid, "%",sep=""), side=3, line=1, cex=1.1)
        dev.off()

    }
}

#########################################################################
