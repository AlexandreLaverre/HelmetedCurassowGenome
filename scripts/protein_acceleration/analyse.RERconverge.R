library(devtools)
library(RERconverge)
library(vioplot)

#######################################################################################
### Defining phenotypes
Helmeted_birds = c("Casuarius_casuarius", "Bucorvus_abyssinicus", "Buceros_rhinoceros",
                   "Pauxi_pauxi", "Numida_meleagris",  "Anseranas_semipalmata", 
                   "Balearica_regulorum", "Anser_cygnoides")
Helmeted_squamates = c("Basiliscus_vittatus", "Chamaeleo_calyptratus")
Helmeted_all_species = c(Helmeted_birds, Helmeted_squamates)

plotRers2 <- function (rermat = NULL, index = NULL, phenv = NULL, rers = NULL, 
                       method = "k", xlims = NULL, plot = 1, xextend = 1, sortrers = F, cex.axis = 10,
                       xlab="Relative Evolution Rate", ylab="Species", colnonIP="deepskyblue3", colIP="brown1") {
  if (is.null(rers)) {
    e1 = rermat[index, ][!is.na(rermat[index, ])]
    colids = !is.na(rermat[index, ])
    e1plot <- e1
    if (exists("speciesNames")) {
      names(e1plot) <- speciesNames[names(e1), ]
    }
    if (is.numeric(index)) {
      gen = rownames(rermat)[index]
    }
    else {
      gen = index
    }
  }
  else {
    e1plot = rers
    gen = "rates"
  }
  names(e1plot)[is.na(names(e1plot))] = ""
  if (!is.null(phenv)) {
    phenvid = phenv[colids]
    fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                       method = method)
    plottitle = paste0(gen, ": rho = ", round(fgdcor$Rho, 3),
                       ", p-val = ", signif(fgdcor$P, digits=3))
    fgd = setdiff(names(e1plot)[phenvid == 1], "")
    df <- data.frame(species = names(e1plot), rer = e1plot, 
                     stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(phenvid > 
                                                                                    0, 2, 1)))
  }
  else {
    plottitle = gen
    fgd = NULL
    df <- data.frame(species = names(e1plot), rer = e1plot, 
                     stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(0, 
                                                                                  2, 1)))
  }
  if (sortrers) {
    df = filter(df, species != "") %>% arrange(desc(rer))
  }
  if (is.null(xlims)) {
    ll = c(min(df$rer) - xextend, max(df$rer) + xextend)
  }
  else {
    ll = xlims
  }
  g <- ggplot(df, aes(x = rer, y = factor(species, levels = unique(ifelse(rep(sortrers, nrow(df)), species[order(rer)], unique(species)))), 
                      col = mole, label = species)) + scale_size_manual(values = c(3,  3)) + geom_point() + 
    scale_color_manual(values = c(colnonIP, colIP)) + scale_x_continuous(limits = ll) + 
    geom_text(hjust = 1.05, size = 3) + ylab(ylab) + xlab(xlab) + 
    ggtitle(plottitle) + geom_vline(xintercept = 0, linetype = "dotted") + 
    
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
          legend.position = "none", panel.background = element_blank(), 
          axis.text = element_text(size = cex.axis, face = "bold", colour = "black"), 
          axis.title = element_text(size = cex.axis, face = "bold"), 
          plot.title = element_text(size = cex.axis, face = "bold")) + 
    theme(axis.line = element_line(colour = "black", size = 1)) + theme(axis.line.y = element_blank())
  
  if (plot) {
    print(g)
  }
  else {
    g
  }
}

#######################################################################################

sp = "all_species" # all_species, birds or squamates
clade = "all"

path = paste("/Users/alaverre/Documents/HelmetedCurassowGenome/results/protein_acceleration", sp, "/", sep="/")

Trees = readRDS(paste0(path, sp, "_RERconverge_Trees.rds"))
RER = readRDS(paste0(path, sp, "_RERconverge_residuals.rds"))
cor = read.table(paste0(path, sp, "_RER_correlations_clade_", clade, ".txt"))

Helmeted = get(paste0("Helmeted_", sp))
pheno <- foreground2Paths(Helmeted, Trees, clade=clade) #terminal, ancestral or all

#######################################################################################
pdf(paste0(path, "RER_results_", sp, ".pdf"), width=6, height=5)
par(mfrow=c(1,2))
par(mai=c(0.5, 0, 0.5, 0))
top_cor = rownames(cor[which(cor$Rho>0),][1,])

##########  Trees ########## 
# Plot Binary tree
binary = foreground2Tree(Helmeted, Trees, clade=clade, plot=F)
bin_tree=plotTreeHighlightBranches(binary, hlspecies=Helmeted, hlcols="red", main="Binary tree")

# Plot top correlated tree
best_tree=plotTreeHighlightBranches(Trees$trees[[top_cor]], hlspecies=Helmeted, hlcols="red", main="Top correlated tree")

# Plot top correlated RER
plotRers2(RER, top_cor, phenv=pheno, sortrers = T, xlims=c(-3, 4))

########## Permutations ########## 
par(mfrow=c(2,2))
par(mai=c(0.7, 0.6, 0.5, 0), mgp=c(2.2,0.8,0))
PermutationsStats <- read.table(paste(path,"/stats_phenotype_permutations.txt",sep=""), h=T)
obsStats <- read.table(paste0(path, "/stats_observed_RER.txt",sep=""), h=T)

# Proportion of positive Rho
obs.val=1-obsStats$Prop.neg
hist(1-PermutationsStats[,"Prop.neg"], xlab="Proportion of positive Rho", main="1000 phenotype permutations", breaks=50, las=1)
abline(v=obs.val, col="red")

# Proportion of pval < 0.05
par(mai=c(0.7, 0.5, 0.5, 0.1))
obs.val = obsStats$nb.signif
hist(PermutationsStats[,"nb.signif"], breaks=40, xlab=paste0("Number of gene (pval<0.05)"), main="", las=1)
abline(v=obs.val, col="red")

##########  Pvalues and FDR per sets ########## 
species = c("all_species", "birds", "squamates")
combined <- c()
N <- c()
for (sp in species){
  path = paste("/Users/alaverre/Documents/HelmetedCurassowGenome/results/protein_acceleration", sp, "/", sep="/")
  cor = read.table(paste0(path, sp, "_RER_correlations_clade_terminal.txt"))
  
  values <- c(cor$P, cor$p.adj)
  type <- c(rep("pval", length(values)/2), rep("FDR", length(values)/2))
  data <- rep(sp, length(values))
  part <- as.data.frame(cbind(values, type, data))
  combined <- rbind(combined, part)
  N <- c(N, paste0("N=",nrow(cor[which(!is.na(cor$P)),])))
  }

combined$values <- as.numeric(combined$values)
combined$type <- factor(combined$type , levels=c("pval", "FDR"))

par(mai=c(0.5, 0.6, 0.5, 0), xpd=T)
vioplot(combined$values ~ combined$type*combined$data, at=c(0,0.8,2.2,3,4.4,5.2),
        col=c("orange2", "mediumseagreen"), border=NA, xlab="", ylab="Values", axes=F, xaxt="n", las=2)
legend("right", c("p-values", "FDR"), fill=c("orange2", "mediumseagreen"), bty='n', inset=-0.35, xpd=NA)
axis(side=1, mgp=c(3, 0.65, 0), at=c(0.4, 2.6, 4.8), labels=c("all species", "birds", "squamates"), cex.axis=1.1)
mtext(N, at=c(0.4, 2.6, 4.8), cex=0.8)

dev.off()

#######################################################################################