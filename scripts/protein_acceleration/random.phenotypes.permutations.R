#!/usr/bin/env Rscript
library(RERconverge)
library(foreach)
library(doSNOW)

#######################################################################################
sp = "all_species" #args[1] # all_species, birds or squamates
clade = "terminal" # considered branches for focal species : all terminal or ancestral

path = paste("/Users/alaverre/Documents/HelmetedCurassowGenome/results/protein_acceleration", sp, "/", sep="/")
Trees = readRDS(paste0(path, sp, "_RERconverge_Trees.rds"))
RER = readRDS(paste0(path, sp, "_RERconverge_residuals.rds"))

### Defining phenotypes
Helmeted_birds = c("Casuarius_casuarius", "Bucorvus_abyssinicus", "Buceros_rhinoceros",
                   "Pauxi_pauxi", "Numida_meleagris",  "Anseranas_semipalmata", 
                   "Balearica_regulorum", "Anser_cygnoides")
Helmeted_squamates = c("Basiliscus_vittatus", "Chamaeleo_calyptratus")
Helmeted_all_species = c(Helmeted_birds, Helmeted_squamates)
Helmeted = get(paste0("Helmeted_", sp))

species <- head(unique(colnames(RER)),-1)
Non.Helmeted = setdiff(species, Helmeted_all_species)

min.sp = ifelse(sp=="squamates", 5, 10)
min.pos = ifelse(sp=="squamates", 2, 5)

##################################################################################

getStatsPermutation <- function(Residuals, Trees, Species){
  phenoIP_random <- foreground2Paths(Species, Trees, clade=clade)
  
  random_cor <- correlateWithBinaryPhenotype(Residuals, phenoIP_random, min.sp=min.sp, min.pos=min.pos, weighted="auto")
  signif <- random_cor[which(random_cor$P < 0.05),]
  
  nb.signif <- nrow(signif)
  med <- median(abs(random_cor$Rho), na.rm=T)
  med.signif <- median(abs(signif$Rho), na.rm=T)
  prop.neg <- nrow(random_cor[which(random_cor$Rho < 0),])/nrow(random_cor)
  prop.neg.signif <- nrow(signif[which(signif$Rho < 0),])/nb.signif
  
  stats <- c(med, med.signif, prop.neg, prop.neg.signif, nb.signif)
  stats
}

##################################################################################
Sys.time()
nb_thread <- 8
nb_permutation <- 1000
cl <- makeCluster(nb_thread)
registerDoSNOW(cl)

pb <- txtProgressBar(max=nb_permutation, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

PermutationsStats <- foreach(i=1:nb_permutation, .combine=rbind, .packages='RERconverge', .options.snow = opts) %dopar% {
  random_species <- sample(species, length(Helmeted))
  tempMatrix = getStatsPermutation(RER, Trees, random_species)
  tempMatrix
}

close(pb)
stopCluster(cl)
Sys.time()

StatsNames = c("Median.Rho", "Median.Rho.signif", "Prop.neg", "Prop.neg.signif", "nb.signif")
colnames(PermutationsStats) <- StatsNames
write.table(PermutationsStats, file=paste0(path, "stats_phenotype_permutations.txt",sep=""), row.names=T, col.names=T, sep="\t",quote=F)

obsStats <- getStatsPermutation(RER, Trees, Helmeted)
write(as.character(StatsNames), file=paste0(path, "stats_observed_RER.txt",sep=""), sep="\t", ncolumns=length(StatsNames))
write(obsStats, file=paste0(path, "stats_observed_RER.txt",sep=""), sep="\t", append=T, ncolumns=length(StatsNames))

##################################################################################
