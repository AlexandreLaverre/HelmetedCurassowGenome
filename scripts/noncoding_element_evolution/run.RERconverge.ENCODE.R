########################################################################

library(RERconverge)

########################################################################

abspath=unlist(strsplit(getwd(), split="\\/"))
lenpath=length(abspath)
path=paste(abspath[1:(lenpath-2)],collapse="/")
pathResults=paste(path,"/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/",sep="")

pathAln=paste(path, "/results/whole_genome_alignments/avian_366/",sep="")

########################################################################

helmeted = c("Casuarius_casuarius", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi", "Numida_meleagris",  "Anseranas_semipalmata",  "Balearica_regulorum", "Anser_cygnoid")

min.sp=10
min.pos=3

########################################################################

for(method in c("iqtree", "phyml")){
  pathTrees=paste(pathResults, method, "_results/",sep="")

  ## prepare input file for RER converge

  if(method=="iqtree"){
    files=system(paste("ls ",pathTrees, " | grep treefile.formatted$",sep=""),intern=T)
  }

  if(method=="phyml"){
    files=system(paste("ls ",pathTrees, " | grep _phyml_tree.formatted$",sep=""),intern=T)
  }

  elements=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))
  names(files)=elements

  est.trees=unlist(lapply(paste(pathTrees,files,sep=""), readLines))
  est.treefiles=data.frame("Gene_name"=as.character(elements), "Newick_tree"=est.trees)

  write.table(est.treefiles, file=paste(pathResults, "RERconverge_",method,"_input_trees.txt",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

########################################################################

  ## read trees

  trees=readTrees(paste(pathResults, "RERconverge_",method,"_input_trees.txt",sep=""))

  ########################################################################

  ## relative rates

  rerw.est=getAllResiduals(trees, transform = "sqrt", weighted=T, scale=T, plot=FALSE)

 ########################################################################

  ## Generating paths for incomplete trees
  pheno <- foreground2Paths(helmeted, trees, clade="all")

  ########################################################################

  ## correlate with phenotype
  cor_all=correlateWithBinaryPhenotype(rerw.est, pheno, min.sp=min.sp, min.pos=min.pos, weighted="auto")

  ordered_cor = cor_all[order(cor_all$p.adj),]

  write.table(ordered_cor, file=paste(pathResults, "RERConverge_correlations_",method,".txt"), row.names=T, col.names=T, sep="\t",quote=F)

}

########################################################################
