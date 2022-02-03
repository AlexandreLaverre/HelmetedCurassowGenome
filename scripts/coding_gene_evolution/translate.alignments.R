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
nbdone=0

for(file in files){
  
  aln=read.fasta(paste(pathResults,"data_for_pelican/",file,sep=""))

  proteins=lapply(aln, translate_alignment)

  write.fasta(proteins, names=names(proteins), file.out=paste(pathResults,"translated_alignments/",file,sep=""))
  
  nbdone=nbdone+1

  if(nbdone%%100==0){
    print(nbdone)
  }
}

########################################################################
  
