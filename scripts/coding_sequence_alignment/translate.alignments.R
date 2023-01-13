########################################################################

library(seqinr)

########################################################################

all.codons=toupper(words(3))
all.aa=unlist(lapply(all.codons, function(x) translate(s2c(x))))
all.aa[which(all.aa=="*")]="X"

all.codons=c(all.codons, "---")
all.aa=c(all.aa, "-")
names(all.aa)=all.codons

translate_alignment <- function(sequence){

  startpos=seq(from=1, to=nchar(sequence), by=3)
  codons=sapply(startpos, function(x) substr(sequence, x, (x+2)))
  protein=paste(all.aa[codons],collapse="")

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
      if(!file.exists(paste(pathResults,"translated_alignments/",file,sep=""))){
        aln=read.fasta(paste(pathResults,"CDS/",file,sep=""), as.string=T, forceDNAtolower=FALSE)

        proteins=lapply(aln, translate_alignment)

        write.fasta(proteins, names=names(proteins), file.out=paste(pathResults,"translated_alignments/",file,sep=""), as.string=T)

        nbdone=nbdone+1

        if(nbdone%%100==0){
          print(nbdone)
        }
      }
    }
  }
}

########################################################################

