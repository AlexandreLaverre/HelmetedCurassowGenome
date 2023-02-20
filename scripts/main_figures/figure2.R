##################################################################

pathFigures="../../results/figures/"

##################################################################

if(load){
    load(paste(pathFigures, "RData/synteny.RData",sep=""))
    load(paste(pathFigures, "RData/gccontent.RData",sep=""))
    load(paste(pathFigures, "RData/repeat.fraction.RData",sep=""))

    load=FALSE
}

##################################################################

if(prepare){

    list.colors=c("black", "red", "indianred", "darkorange", "steelblue", "blue", "seagreen", "darkolivegreen", "purple", "lightpink", "yellow", "slategray2", "turquoise")

    selected.chr=list()
    selected.synteny=list()
    ref.colors=list()
    legend.chr=list()
    ref.sp=c("Gallus_gallus", "Anolis_carolinensis")
    names(ref.sp)=c("Pauxi_pauxi", "Basiliscus_vittatus")

    for(sp in c("Basiliscus_vittatus", "Pauxi_pauxi")){
        this.gcinfo=gcinfo[[sp]]
        this.gcinfo=this.gcinfo[order(this.gcinfo$Size, decreasing=T),]
        selected.chr[[sp]]=this.gcinfo[which(this.gcinfo$Size>=20e6),]

        this.synteny=synteny[[sp]]
        this.synteny=this.synteny[which(this.synteny$Chr%in%selected.chr[[sp]][,"Chr"]),]
        this.synteny$MeanPos=(this.synteny$Start+this.synteny$End)/2


        nb.ref.chr=as.numeric(table(this.synteny$ReferenceChr))
        names(nb.ref.chr)=levels(as.factor(this.synteny$ReferenceChr))
        selected.ref.chr=names(nb.ref.chr)[which(nb.ref.chr>20)]

        print(length(selected.ref.chr))

        selected.colors=list.colors[1:length(selected.ref.chr)]
        names(selected.colors)=selected.ref.chr

        autoassembled=intersect(selected.ref.chr, as.character(1:50))
        other=setdiff(selected.ref.chr,autoassembled)
        order.chr=c(autoassembled[order(as.numeric(autoassembled))],other)
        legend.chr[[sp]]=selected.colors[order.chr]

        this.synteny$Color="gray60"
        this.synteny$Color[which(this.synteny$ReferenceChr%in%selected.ref.chr)]=selected.colors[this.synteny$ReferenceChr[which(this.synteny$ReferenceChr%in%selected.ref.chr)]]

        selected.synteny[[sp]]=this.synteny
    }

    prepare=FALSE
}

##################################################################

## page width 8.5, page height 11

pdf(paste(pathFigures, "main_figures/Figure2.pdf",sep=""), width=7.5, height=6.5)

##################################################################

## figure layout

m=matrix(rep(NA, 10*10),nrow=10)

for(i in 1:6){
    m[i,]=c(rep(1,5), rep(2,5))
}

for(i in 7:10){
    m[i,]=c(rep(3,10))
}

layout(m)

##################################################################

## first plots: synteny with other species

labels=c("A", "B")
names(labels)=c("Pauxi_pauxi", "Basiliscus_vittatus")

for(sp in c("Pauxi_pauxi", "Basiliscus_vittatus")){

    this.ref=paste(unlist(strsplit(ref.sp[sp], split="_")), collapse=" ")

    this.synteny=selected.synteny[[sp]]

    this.chr=selected.chr[[sp]]
    nb.chr=nrow(this.chr)
    max.size=max(this.chr$Size)

    smally=max.size/20
    ylim=c(-smally, max.size+smally)
    xlim=c(0.5, nb.chr+0.5)

    par(mar=c(2.1, 3.75, 3.1, 1.1))
    plot(1, type="n", axes=F, xlab="", ylab="", xlim=xlim, ylim=ylim)

    for(i in 1:nb.chr){
        this.chr.name=this.chr$Chr[i]
        segments(i, 1, i, this.chr$Size[i], col="gray40")

        this.synteny.chr=this.synteny[which(this.synteny$Chr==this.chr.name),]

        other=which(this.synteny.chr$Color=="gray60")
        sel=which(this.synteny.chr$Color!="gray60")

        points(rep(i, length(other)), this.synteny.chr$MeanPos[other], col=this.synteny.chr$Color[other], pch=20)
        points(rep(i, length(sel)), this.synteny.chr$MeanPos[sel], col=this.synteny.chr$Color[sel], pch=20)

        short.name=unlist(strsplit(this.chr.name, split="Scaffold"))[2]

        text(short.name, x=i, y=-smally)
    }

    this.legend=legend.chr[[sp]]

    legend("topright", legend=names(this.legend), col=this.legend, pch=20, bty="n", inset=c(0.05,0.05))

    if(sp=="Pauxi_pauxi"){
        text(this.ref, font=3, x=10, y=max.size*2/3, srt=90, cex=1.1)
    }

    if(sp=="Basiliscus_vittatus"){
        text(this.ref, font=3, x=8, y=max.size*2.35/3, srt=90, cex=1.1)
    }

    yax=pretty(c(1, max.size/1e6))
    yax=yax[which(yax<(ylim[2]/1e6))]
    yaxlab=paste(yax, "Mb",sep="")
    axis(side=2, at=yax*1e6, labels=yaxlab, mgp=c(3,0.5,0))

    synsp=paste(unlist(strsplit(sp, split="_")), collapse=" ")
    mtext(synsp, font=3, side=3, line=1, cex=0.8)

    mtext(labels[sp], font=2, side=3, at=-1.5, line=0.5, cex=1.1)

    mtext("scaffold size", side=2, line=2.5, cex=0.75)
}

##################################################################

dev.off()

##################################################################
