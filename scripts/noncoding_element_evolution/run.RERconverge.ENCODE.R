########################################################################

library(RERconverge)

########################################################################

abspath=unlist(strsplit(getwd(), split="\\/"))
lenpath=length(abspath)
path=paste(abspath[1:(lenpath-2)],collapse="/")
pathResults=paste(path,"/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/",sep="")
pathTrees=paste(pathResults, "iqtree_results/",sep="")
pathAln=paste(path, "/results/whole_genome_alignments/avian_366/",sep="")

########################################################################

## prepare input file for RER converge

files=system(paste("ls ",pathTrees, " | grep .treefile$",sep=""),intern=T)

elements=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\.treefile"))[1]))

trees=unlist(lapply(paste(pathTrees,files,sep=""), readLines))

treefiles=data.frame("Gene_name"=as.character(elements), "Newick_tree"=trees)

write.table(treefiles[1:100,], file=paste(pathResults, "RERconverge_trees.txt",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

########################################################################

masterTree=read.tree(paste(pathAln, "Birds366_tree.txt",sep=""))

spList=c("Struthio_camelus", "Dromaius_novaehollandiae", "Gallus_gallus", "Meleagris_gallopavo", "Penelope_pileata", "Alectura_lathami", "Anas_platyrhynchos_platyrhynchos", "Grus_americana", "Calidris_pugnax", "Upupa_epops", "Rhinopomastus_cyanomelas", "Strix_occidentalis", "Ficedula_albicollis", "Aquila_chrysaetos", "Geospiza_fortis", "Parus_major", "Serinus_canaria", "Coturnix_japonica", "Casuarius_casuarius", "Numida_meleagris", "Anseranas_semipalmata", "Anser_cygnoid", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

masterTree=keep.tip(masterTree, spList)

########################################################################

## readtrees

trees=readTrees(paste(pathResults, "RERconverge_trees.txt",sep=""), masterTree=masterTree, reestimateBranches=TRUE)

########################################################################

## relative rates

rerw=getAllResiduals(trees, transform = "sqrt", weighted = T, scale = T, plot=F)

########################################################################
