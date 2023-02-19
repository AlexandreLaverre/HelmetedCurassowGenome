################################################################

pathOrtho="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/Orthologues/"
pathAnnot="../../results/genome_annotation/"
pathEnsembl="../../data/ensembl_annotations/"
pathResults="../../results/synteny_analysis/"

release=103

################################################################

species=c("Pauxi_pauxi", "Basiliscus_vittatus")
refspecies.annot=c("Chicken", "Anolis_carolinensis")
refspecies.ortho=c("Gallus_gallus", "Anolis_carolinensis")

################################################################

for(i in 1:length(species)){
    this.sp=species[i]
    this.refsp.annot=refspecies.annot[i]
    this.refsp.ortho=refspecies.ortho[i]

    ## read orthologues

    ortho=read.table(paste(pathOrtho, "Orthologues_",this.sp,"/",this.sp,"__v__",this.refsp.ortho,".tsv", sep=""),h=T, stringsAsFactors=F, sep="\t")

    ## keep only one to one
    multi=apply(ortho[,c(this.sp, this.refsp.ortho)],1, function(x) length(grep(",",x))>0)
    ortho=ortho[!multi,c(this.sp, this.refsp.ortho)]

    ## change gene id to match gene coordinates
    ortho[,this.refsp.ortho]=unlist(lapply(ortho[,this.refsp.ortho], function(x) unlist(strsplit(x, split="\\."))[1]))

    ## gene coords in ref sp
    coords.ref=read.table(paste(pathEnsembl, this.refsp.annot, "/GeneCoordinates_Ensembl", release, ".txt",sep=""),h=T, stringsAsFactors=F, sep="\t", quote="\"")
    colnames(coords.ref)=c("GeneID", "TranscriptID", "ProteinID", "Chr", "Start", "End", "Strand")
    coords.ref=coords.ref[which(!duplicated(coords.ref$GeneID)), c("GeneID", "Chr", "Start", "End", "Strand")]
    rownames(coords.ref)=coords.ref$GeneID

    ## gene coords in target sp
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

    ## assemble all info for ortho pairs
    ortho.info=data.frame("GeneID"=ortho[,this.sp], "Chr"=coords.sp[ortho[,this.sp], "Chr"], "Start"=coords.sp[ortho[,this.sp], "Start"], "End"=coords.sp[ortho[,this.sp], "End"], "Strand"=coords.sp[ortho[,this.sp], "Strand"], "ReferenceGeneID"=ortho[,this.refsp.ortho], "ReferenceChr"=coords.ref[ortho[,this.refsp.ortho], "Chr"], "ReferenceStart"=coords.ref[ortho[,this.refsp.ortho], "Start"], "ReferenceEnd"=coords.ref[ortho[,this.refsp.ortho], "End"], "ReferenceStrand"=coords.ref[ortho[,this.refsp.ortho], "Strand"])

    ortho.info=ortho.info[order(ortho.info$Start),]
    ortho.info=ortho.info[order(ortho.info$Chr),]

    if(dir.exists(paste(pathResults, this.sp, sep=""))){
        print("output dir exists")
    } else{
        system(paste("mkdir -p ",pathResults, this.sp, sep=""))
    }

    write.table(ortho.info, file=paste(pathResults, this.sp, "/OrthoGeneCoordinates_",this.refsp.ortho,".txt", sep=""), row.names=F, col.names=T, sep="\t")
}

################################################################
