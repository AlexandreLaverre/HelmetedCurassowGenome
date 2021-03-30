####################################################################

pathAnnot="../../data/ensembl_annotations/"
pathResults="../../results/genome_assembly_quality/MEGAHIT_RAGOUT/"

####################################################################

ensrelease=103
minpcid=60

####################################################################

for(ref in c("Chicken", "Duck")){

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

}

####################################################################
