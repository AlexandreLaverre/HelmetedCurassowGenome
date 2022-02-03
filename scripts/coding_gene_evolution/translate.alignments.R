########################################################################

library(seqinr)

########################################################################

pathResults="../../results/coding_gene_evolution/"

########################################################################

files=system(paste("ls ",pathResults, "data_for_pelican/", sep=""), intern=T)

########################################################################

translate_alignment <- function(seq){
  startpos=seq(from=1, to=length(seq), by=3)

  protein=unlist(lapply(startpos, function(x) {codon=seq[x:(x+2)]; if(all(codon=="-")){return("-")} else{return(translate(codon))}}))

  return(protein)
}

########################################################################

results=read.table(paste(pathResults, "pelican_output_by_category/best_sites.tsv",sep=""),h=T, stringsAsFactors=F)

########################################################################

nbdone=0

for(gene in unique(results$alignment)){
  file=paste(gene, ".fa", sep="")
  
  aln=read.fasta(paste(pathResults,"data_for_pelican/",file,sep=""))

  proteins=lapply(aln, translate_alignment)

  write.fasta(proteins, names=names(proteins), file.out=paste(pathResults,"translated_alignments/",file,sep=""))
  
  nbdone=nbdone+1

  if(nbdone%%100==0){
    print(nbdone)
  }
}

########################################################################
  
