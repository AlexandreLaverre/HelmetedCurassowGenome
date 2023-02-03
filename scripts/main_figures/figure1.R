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

    for(sp in helmeted.species){
        sp.images[[sp]]=readPNG(paste(pathImages,sp,"_head_color.png",sep=""))
    }

    sp.pos=list()
    sp.pos[["Numida_meleagris"]]=list("xleft"=0.35, "ybottom"=35.5, "xright"=0.42, "ytop"=39,xpd=NA)
    sp.pos[["Pauxi_pauxi"]]=list("xleft"=0.45, "ybottom"=32, "xright"=0.52, "ytop"=35,xpd=NA)
    sp.pos[["Anser_cygnoides"]]=list("xleft"=0.36, "ybottom"=29, "xright"=0.43, "ytop"=33,xpd=NA)
    sp.pos[["Anseranas_semipalmata"]]=list("xleft"=0.43, "ybottom"=25, "xright"=0.54, "ytop"=29,xpd=NA)
    sp.pos[["Buceros_rhinoceros"]]=list("xleft"=0.36, "ybottom"=21, "xright"=0.45, "ytop"=25,xpd=NA)
    sp.pos[["Bucorvus_abyssinicus"]]=list("xleft"=0.45, "ybottom"=18, "xright"=0.55, "ytop"=22.5,xpd=NA)
    sp.pos[["Balearica_regulorum"]]=list("xleft"=0.34, "ybottom"=14, "xright"=0.44, "ytop"=19,xpd=NA)
    sp.pos[["Casuarius_casuarius"]]=list("xleft"=0.47, "ybottom"=12, "xright"=0.54, "ytop"=17,xpd=NA)
    sp.pos[["Basiliscus_vittatus"]]=list("xleft"=0.45, "ybottom"=7.25, "xright"=0.55, "ytop"=11,xpd=NA)
    sp.pos[["Chamaeleo_calyptratus"]]=list("xleft"=0.46, "ybottom"=1.5, "xright"=0.55, "ytop"=6,xpd=NA)

    label.pos=list()
    label.pos[["Numida_meleagris"]]=0.29
    label.pos[["Pauxi_pauxi"]]=0.25
    label.pos[["Anser_cygnoides"]]=0.27
    label.pos[["Anseranas_semipalmata"]]=0.305
    label.pos[["Buceros_rhinoceros"]]=0.33
    label.pos[["Bucorvus_abyssinicus"]]=0.305
    label.pos[["Balearica_regulorum"]]=0.305
    label.pos[["Casuarius_casuarius"]]=0.275
    label.pos[["Basiliscus_vittatus"]]=0.385
    label.pos[["Chamaeleo_calyptratus"]]=0.395

    yseg=list()
    yseg[["Numida_meleagris"]]=37.5
    yseg[["Pauxi_pauxi"]]=33.5
    yseg[["Anser_cygnoides"]]=31
    yseg[["Anseranas_semipalmata"]]=26.5
    yseg[["Buceros_rhinoceros"]]=22.5
    yseg[["Bucorvus_abyssinicus"]]=19
    yseg[["Balearica_regulorum"]]=16
    yseg[["Casuarius_casuarius"]]=14
    yseg[["Basiliscus_vittatus"]]=9
    yseg[["Chamaeleo_calyptratus"]]=5.25


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

pdf(paste(pathFigures, "main_figures/Figure1.pdf",sep=""), width=6.5, height=6.5)

## figure layout

m=matrix(rep(NA, 10*10),nrow=10)

for(i in 1:10){
    m[i,]=rep(1,10)
}

layout(m)

##################################################################

## first plot: species phylogeny

par(mar=c(0.5, 0.5, 0.5, 0.05))

plot(sptree, x.lim=c(0,0.65), label.offset=0.0025, edge.color=edge.color, edge.width=1.5, font=tip.font, cex=1.1)

add.scale.bar(x=0.01, y=0.5)

##################################################################

## second plot: species images

xshift=0.1
xshiftlab=0.035

for(sp in names(sp.pos)){

    rasterImage(sp.images[[sp]], "xleft"=xshift+sp.pos[[sp]][["xleft"]], "xright"=xshift+sp.pos[[sp]][["xright"]], "ybottom"=sp.pos[[sp]][["ybottom"]], "ytop"=sp.pos[[sp]][["ytop"]])

    tip.index=which(sptree$tip.label==sp)

    if(sp=="Anser_cygnoides"){
        tip.index=tip.index+1
    }

    if(sp=="Chamaeleo_calyptratus"){
        tip.index=tip.index-1
    }

    this.yseg=yseg[[sp]]

    xseg=sp.pos[[sp]][["xleft"]]+0.01

    if(sp%in%c("Pauxi_pauxi", "Numida_meleagris")){
        xseg=sp.pos[[sp]][["xleft"]]-0.01
    }

    if(sp%in%c("Casuarius_casuarius")){
        xseg=sp.pos[[sp]][["xleft"]]-0.02
    }
     if(sp%in%c("Anseranas_semipalmata", "Basiliscus_vittatus")){
        xseg=sp.pos[[sp]][["xleft"]]+0.02
     }

    if(sp%in%c("Chamaeleo_calyptratus")){
        xseg=sp.pos[[sp]][["xleft"]]+0.04
    }


    segments(x0=label.pos[[sp]]+xshiftlab, x1=xshift+xseg, y0=tip.index, y1=this.yseg, lty=3, xpd=NA)
}

##################################################################

dev.off()

##################################################################
