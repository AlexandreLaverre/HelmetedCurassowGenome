########################################################################

library(RERconverge)

method="iqtree"

########################################################################

abspath=unlist(strsplit(getwd(), split="\\/"))
lenpath=length(abspath)
path=paste(abspath[1:(lenpath-2)],collapse="/")
pathResults=paste(path,"/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/",sep="")
pathTrees=paste(pathResults, method, "_results/",sep="")
pathAln=paste(path, "/results/whole_genome_alignments/avian_366/",sep="")

########################################################################

masterTree=read.tree(paste(pathAln, "Birds366_tree.txt",sep=""))

spList=c("Struthio_camelus", "Dromaius_novaehollandiae", "Gallus_gallus", "Meleagris_gallopavo", "Penelope_pileata", "Alectura_lathami", "Anas_platyrhynchos_platyrhynchos", "Grus_americana", "Calidris_pugnax", "Upupa_epops", "Rhinopomastus_cyanomelas", "Strix_occidentalis", "Ficedula_albicollis", "Aquila_chrysaetos", "Geospiza_fortis", "Parus_major", "Serinus_canaria", "Coturnix_japonica", "Casuarius_casuarius", "Numida_meleagris", "Anseranas_semipalmata", "Anser_cygnoid", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

masterTree=keep.tip(masterTree, spList)

masterTree$node.label=NULL

########################################################################

## prepare input file for RER converge

if(method=="iqtree"){
    files=system(paste("ls ",pathTrees, " | grep treefile.formatted$",sep=""),intern=T)
}

if(method=="phyml"){
    files=system(paste("ls ",pathTrees, " | grep _phyml_tree.formatted$",sep=""),intern=T)
}

elements=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))
names(files)=elements

for(el in elements){
    this.tree=read.tree(file=paste(pathTrees,files[el],sep=""))
    this.tree$node.label=NULL
}

est.trees=unlist(lapply(paste(pathTrees,files,sep=""), readLines))
est.treefiles=data.frame("Gene_name"=as.character(elements), "Newick_tree"=est.trees)

write.table(est.treefiles, file=paste(pathResults, "RERconverge_",method,"_input_trees.txt",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

########################################################################

## readtrees

trees=readTrees(paste(pathResults, "RERconverge_",method,"_input_trees.txt",sep=""))

########################################################################

## relative rates

rerw.est=getAllResiduals(est.trees, transform = "sqrt", weighted=T, scale=T)

########################################################################
########################################################################



########################################################################
########################################################################
