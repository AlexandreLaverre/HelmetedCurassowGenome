this.tree=multi2di(this.tree)

write.tree(this.tree,paste(pathTrees, files[el],".formatted",sep=""))


tr1=read.tree(paste(pathTrees,files[1],".formatted",sep=""))
tr1$node.label=NULL

tr2=keep.tip(masterTree, tr1$tip.label)
tr2$node.label=NULL

tr1=root(tr1, outgroup="Gallus_gallus")
tr2=root(tr2, outgroup="Gallus_gallus")


newtree1=tr2

nbsp=length(tr1$tip.label)

edge1=tr1$edge
edge2=tr2$edge

desc1nodes=unlist(lapply(unique(edge1[,1]), function(x) paste(sort(extract.clade(tr1, x)$tip.label), collapse=",")))
desc2nodes=unlist(lapply(unique(edge2[,1]), function(x) paste(sort(extract.clade(tr2, x)$tip.label), collapse=",")))

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

print(all(br1%in%br2))

len1=tr1$edge.length

names(len1)=br1

newtree1$edge.length=as.numeric(len1[br2])
