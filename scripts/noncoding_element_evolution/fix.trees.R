###################################################################################

library(ape)

###################################################################################

abspath=unlist(strsplit(getwd(), split="\\/"))
lenpath=length(abspath)
path=paste(abspath[1:(lenpath-2)],collapse="/")
pathResults=paste(path,"/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/",sep="")
pathAln=paste(path, "/results/whole_genome_alignments/avian_366/",sep="")

########################################################################

masterTree=read.tree(paste(pathAln, "Birds366_tree.txt",sep=""))

spList=c("Struthio_camelus", "Dromaius_novaehollandiae", "Gallus_gallus", "Meleagris_gallopavo", "Penelope_pileata", "Alectura_lathami", "Anas_platyrhynchos_platyrhynchos", "Grus_americana", "Calidris_pugnax", "Upupa_epops", "Rhinopomastus_cyanomelas", "Strix_occidentalis", "Ficedula_albicollis", "Aquila_chrysaetos", "Geospiza_fortis", "Parus_major", "Serinus_canaria", "Coturnix_japonica", "Casuarius_casuarius", "Numida_meleagris", "Anseranas_semipalmata", "Anser_cygnoid", "Balearica_regulorum", "Bucorvus_abyssinicus", "Buceros_rhinoceros", "Pauxi_pauxi")

masterTree=keep.tip(masterTree, spList)

########################################################################

get.node.desc<-function(tree,x){
  cl=extract.clade(tree, x)
  descsp1=paste(sort(cl$tip.label), collapse=",")
  descsp2=paste(sort(setdiff(tree$tip.label,cl)), collapse=",")
  id=paste(sort(c(descsp1, descsp2)), collapse="-")
  return(id)
}

########################################################################

for(method in c("iqtree", "phyml")){

  pathTrees=paste(pathResults, method, "_results/",sep="")

  if(method=="iqtree"){
    files=system(paste("ls ",pathTrees, " | grep treefile$",sep=""),intern=T)
  }

  if(method=="phyml"){
    files=system(paste("ls ",pathTrees, " | grep _phyml_tree.txt",sep=""),intern=T)
  }

  elements=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))
  names(files)=elements

  for(el in elements){
    print(el)

    tr1=read.tree(file=paste(pathTrees,files[el],sep=""))
    tr1$node.label=NULL
    tr1=multi2di(tr1)
    tr1=root(tr1, outgroup="Gallus_gallus")

    tr2=keep.tip(masterTree, tr1$tip.label)
    tr2$node.label=NULL
    tr2=root(tr2, outgroup="Gallus_gallus")

    if(as.numeric(dist.topo(tr1, tr2)==0)){

      newtree1=tr2

      nbsp=length(tr1$tip.label)

      edge1=tr1$edge
      edge2=tr2$edge

      desc1nodes=unlist(lapply(unique(edge1[,1]), function(x) get.node.desc(tr1,x)))
      desc2nodes=unlist(lapply(unique(edge2[,1]), function(x) get.node.desc(tr2,x)))

      names(desc1nodes)=as.character(unique(edge1[,1]))
      names(desc2nodes)=as.character(unique(edge2[,1]))

      desc1tip=tr1$tip.label
      names(desc1tip)=as.character(1:nbsp)

      desc2tip=tr2$tip.label
      names(desc2tip)=as.character(1:nbsp)

      desc1=c(desc1nodes, desc1tip)
      desc2=c(desc2nodes, desc2tip)

      br1=paste(desc1[as.character(edge1[,1])], "-",desc1[as.character(edge1[,2])])
      br2=paste(desc2[as.character(edge2[,1])], "-",desc2[as.character(edge2[,2])])

      if(all(br1%in%br2) & all(br2%in%br1)){

        len1=tr1$edge.length

        names(len1)=br1

        newtree1$edge.length=as.numeric(len1[br2])

        write.tree(newtree1, file=paste(pathTrees, files[el], ".formatted", sep=""))
      } else{
        stop("weird, branch lengths do not correspond")
      }
    } else{
      print("weird, tree topologies do not correspond")
    }
  }
}

###################################################################################
