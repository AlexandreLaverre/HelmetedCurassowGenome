########################################################################

library(seqinr)

########################################################################

translate_alignment <- function(seq){
  startpos=seq(from=1, to=length(seq), by=3)

  protein=unlist(lapply(startpos, function(x) {codon=seq[x:(x+2)]; if(all(codon=="-")){return("-")} else{return(translate(codon))}}))

  return(protein)
}

########################################################################

for(spset in c("all_species", "without_chameleons")){
  for(dataset in c("all_species", "birds", "squamates")){
    
    pathResults=paste("../../results/coding_gene_evolution/",spset, "/", dataset, "/",sep="")

    if(dir.exists(paste(pathResults,"translated_alignments/",sep=""))){
      print("path output already there")
    } else{
      system(paste("mkdir ",pathResults,"translated_alignments/",sep=""))
    }
    
    ########################################################################
    
    files=system(paste("ls ",pathResults, "CDS/ | grep aln.best.fas", sep=""), intern=T)
    
    ########################################################################
    
    nbdone=0
    
    for(file in files){
            
      aln=read.fasta(paste(pathResults,"CDS/",file,sep=""))
      
      proteins=lapply(aln, translate_alignment)
      
      write.fasta(proteins, names=names(proteins), file.out=paste(pathResults,"translated_alignments/",file,sep=""))
      
      nbdone=nbdone+1
      
      if(nbdone%%100==0){
        print(nbdone)
      }
    }
  }
}

########################################################################
  
