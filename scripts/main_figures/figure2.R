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

    ## proportion repeats


    allspecies=c("Pauxi_pauxi", "Basiliscus_vittatus", "Gallus_gallus", "Anolis_carolinensis")

    pc.rep=list()

    for(sp in allspecies){
        this.gcinfo=gcinfo[[sp]]
        tot.size=sum(this.gcinfo$Size)

        this.rep=repfr[[sp]]

        selclasses=list("DNA", "SINE", "LINE", "LTR", c("Low_complexity", "Simple_repeat"))
        selclasses[[length(selclasses)+1]]=setdiff(this.rep$RepeatClass, unlist(selclasses))
        names(selclasses)=c("DNA", "SINE", "LINE", "LTR", "Low/Simple", "Other")

        this.pc.rep=unlist(lapply(selclasses, function(x) sum(this.rep$TotalLength[which(this.rep$RepeatClass%in%x)])))
        this.pc.rep=100*this.pc.rep/tot.size

        pc.rep[[sp]]=this.pc.rep
    }

    prepare=FALSE
}

##################################################################

## page width 8.5, page height 11

pdf(paste(pathFigures, "main_figures/Figure2.pdf",sep=""), width=7.5, height=6.5)

##################################################################

## figure layout

m=matrix(rep(NA, 10*12),nrow=10)

for(i in 1:6){
    m[i,]=c(rep(1, 6), rep(2, 6))
}

for(i in 7:10){
    m[i,]=c(rep(3, 4), rep(5, 2), rep(4, 4), rep(6, 2))
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
    xlim=c(0.75, nb.chr+0.25)

    par(mar=c(2.1, 3.75, 2.5, 1.1))
    plot(1, type="n", axes=F, xlab="", ylab="", xlim=xlim, ylim=ylim)

    for(i in 1:nb.chr){
        this.chr.name=this.chr$Chr[i]
        segments(i, 1, i, this.chr$Size[i], col="gray40")

        this.synteny.chr=this.synteny[which(this.synteny$Chr==this.chr.name),]

        other=which(this.synteny.chr$Color=="gray60")
        sel=which(this.synteny.chr$Color!="gray60")

        points(rep(i, length(other)), this.synteny.chr$MeanPos[other], col=this.synteny.chr$Color[other], pch=20, xpd=NA)
        points(rep(i, length(sel)), this.synteny.chr$MeanPos[sel], col=this.synteny.chr$Color[sel], pch=20, xpd=NA)

        short.name=unlist(strsplit(this.chr.name, split="Scaffold"))[2]

        text(short.name, x=i, y=-smally, cex=1.05)
    }

    this.legend=legend.chr[[sp]]

    legend("topright", legend=names(this.legend), col=this.legend, pch=20, bty="n", inset=c(0.05,0.05))

    if(sp=="Pauxi_pauxi"){
        text(this.ref, font=3, x=10, y=max.size*2.1/3, srt=90, cex=1.1)
    }

    if(sp=="Basiliscus_vittatus"){
        text(this.ref, font=3, x=8, y=max.size*2.455/3, srt=90, cex=1.1)
    }

    yax=pretty(c(1, max.size/1e6))
    yax=yax[which(yax<(ylim[2]/1e6))]
    yaxlab=paste(yax, "Mb",sep="")
    axis(side=2, at=yax*1e6, labels=yaxlab, mgp=c(3,0.5,0))

    synsp=paste(unlist(strsplit(sp, split="_")), collapse=" ")
    mtext(synsp, font=3, side=3, line=1, cex=0.8)

    mtext(labels[sp], font=2, side=3, at=-1.25, line=1, cex=1.1)

    mtext("scaffold size", side=2, line=2.5, cex=0.75)
}

##################################################################

## second and third plot: GC content as a function of chromosome size

labels=c("C", "E")
names(labels)=c("Pauxi_pauxi", "Basiliscus_vittatus")

for(sp in c("Pauxi_pauxi", "Basiliscus_vittatus")){
    this.gcinfo=gcinfo[[sp]]
    this.gcinfo=this.gcinfo[which(this.gcinfo$Size>5e6),]

    other.gcinfo=gcinfo[[ref.sp[[sp]]]]
    other.gcinfo=other.gcinfo[which(other.gcinfo$Size>5e6),]

    par(mar=c(3.5,3.5,2.1,0.5))

    ylim=range(c(this.gcinfo$GC, other.gcinfo$GC))
    smally=5*diff(ylim)/100
    ylim=ylim+c(-smally, +smally)
    xlim=range(c(this.gcinfo$Size, other.gcinfo$GC))+c(0, 30e6)

    plot(this.gcinfo$Size, this.gcinfo$GC, pch=20, xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, col="red")
    points(other.gcinfo$Size, other.gcinfo$GC, col="black", pch=20)

    xax=pretty(this.gcinfo$Size/1e6)
    xax=xax[which((xax*1e6)<=xlim[2])]
    axis(side=1, at=xax*1e6, labels=paste(xax), mgp=c(3,0.5,0))
    axis(side=2,  mgp=c(3,0.65,0))

    mtext("scaffold size (Mb)", side=1, line=2, cex=0.75)
    mtext("scaffold GC content", side=2, line=2.1, cex=0.75)


    mtext(labels[sp], font=2, side=3, at=-55e6, line=1, cex=1.1)


    this.ref=paste(unlist(strsplit(ref.sp[sp], split="_")), collapse=" ")
    synsp=paste(unlist(strsplit(sp, split="_")), collapse=" ")

    legend("topright", legend=c(synsp, this.ref), pch=20, col=c("red", "black"), bty="n", inset=0.01, text.font=3)
}

##################################################################

## 4th and 5th plot: repeat type

labels=c("D", "F")
names(labels)=c("Pauxi_pauxi", "Basiliscus_vittatus")


##################################################################

dev.off()

##################################################################
