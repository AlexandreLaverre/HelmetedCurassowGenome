#########################################################################

library(ape)
library(seqinr)

#########################################################################

args = commandArgs(trailingOnly=TRUE)

#########################################################################

path.species.tree=args[1]
path.aln=args[2]
path.output.tree=args[3]

#########################################################################

sp.tree=read.tree(file=path.species.tree)
sp.tree$node.label <- NULL

aln=read.fasta(path.aln)

this.tree=keep.tip(phy=sp.tree, names(aln))
write.tree(this.tree, file=path.output.tree)

#########################################################################

