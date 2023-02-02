##################################################################

library(ape)
library(png)

##################################################################

pathFigures="../../results/figures/"
pathImages="../../docs/species_images/"

##################################################################

if(load){
    load(paste(pathFigures, "RData/species.phylogeny.RData",sep=""))

    load=FALSE
}

##################################################################

if(prepare){

    ## helmeted species

    helmeted.species=c("Numida_meleagris", "Pauxi_pauxi", "Anser_cygnoides", "Anseranas_semipalmata", "Buceros_rhinoceros", "Bucorvus_abyssinicus", "Balearica_regulorum", "Casuarius_casuarius", "Basiliscus_vittatus", "Chamaeleo_calyptratus")

    ## species images
    sp.images=list()
    sp.pos=list()

    for(sp in helmeted.species){
        sp.images[[sp]]=readPNG(paste(pathImages,sp,".png",sep=""))
    }

    sp.pos[["Numida_meleagris"]]=list("xleft"=0.1, "ybottom"=37, "xright"=0.75, "ytop"=38)

    ## rotate phylogeny to accomodate species images

    node1=getMRCA(sptree, c("Anser_cygnoides", "Anser_brachyrhynchus", "Anas_platyrhynchos_platyrhynchos"))
    sptree=rotate(sptree, node=node1)


    node2=getMRCA(sptree, c("Pogona_vitticeps", "Chamaeleo_calyptratus", "Chamaeleo_chamaeleon_recticrista"))
    sptree=rotate(sptree, node=node2)

    index.tips=which(sptree$tip.label%in%helmeted.species)

    edge.color=rep("gray40", nrow(sptree$edge))
    edge.color[which(sptree$edge[,2]%in%index.tips)]="red"

    index.bucanc=getMRCA(sptree, c("Buceros_rhinoceros", "Bucorvus_abyssinicus"))
    edge.color[which(sptree$edge[,2]%in%index.bucanc)]="red"

    tip.font=rep(3, length(sptree$tip.label))
    tip.font[which(sptree$tip.label%in%c("Pauxi_pauxi", "Basiliscus_vittatus", "Chamaeleo_calyptratus", "Chamaeleo_chamaeleon_recticrista"))]=4


    prepare=FALSE
}

##################################################################

## page width 8.5, page height 11

pdf(paste(pathFigures, "main_figures/Figure1.pdf",sep=""), width=5.1, height=6.5)

## figure layout

m=matrix(rep(NA, 10*10),nrow=10)

for(i in 1:10){
    m[i,]=c(rep(1,7), rep(2, 3))
}

layout(m)

##################################################################

## first plot: species phylogeny

par(mar=c(0.5, 0.5, 0.5, 0))

plot(sptree, x.lim=c(0,0.7), label.offset=0.0025, edge.color=edge.color, edge.width=1.5, font=tip.font, cex=1.1)

##################################################################

## second plot: species images


par(mar=c(0.5, 0, 0.5, 0.5))

plot(1, type="n", axes=F,  xlim=c(0,1), ylim=c(0,1))


nbsp=length(sptree$tip.label)

for(sp in c("Numida_meleagris")){
    nbpix.x=dim(sp.images[[sp]])[1]
    nbpix.y=dim(sp.images[[sp]])[2]

    height=1.5
    width=nbpix.x/nbpix.y/height

    sp.pos[[sp]]=list("xleft"=0.2, "ybottom"=(nbsp-3)/nbsp, "xright"=0.2+width, "ytop"=1.03,xpd=NA)

    rect("xleft"=sp.pos[[sp]][["xleft"]], "xright"=sp.pos[[sp]][["xright"]], "ybottom"=sp.pos[[sp]][["ybottom"]], "ytop"=sp.pos[[sp]][["ytop"]])


    rasterImage(sp.images[[sp]], "xleft"=sp.pos[[sp]][["xleft"]], "xright"=sp.pos[[sp]][["xright"]], "ybottom"=sp.pos[[sp]][["ybottom"]], "ytop"=sp.pos[[sp]][["ytop"]])
}

##################################################################

dev.off()

##################################################################
