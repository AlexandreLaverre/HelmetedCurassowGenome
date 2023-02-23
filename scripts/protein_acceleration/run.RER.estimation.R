#!/usr/bin/env Rscript

#######################################################################################
# Load and install packages if required
list.of.packages <- c("devtools", "RColorBrewer", "gplots", "phytools", "ape", "maps", "Rcpp","geiger", "knitr", "RcppArmadillo",  "weights", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("ggtree")) BiocManager::install("ggtree")

library(devtools)
if (!require("RERconverge")) install_github("nclark-lab/RERconverge")
try(library(RERconverge), silent=TRUE)

#######################################################################################
# Define path
#args = commandArgs(trailingOnly=TRUE)

sp = "birds" #args[1] # all_species, birds or squamates
clade = "all" # considered branches for focal species : all terminal or ancestral

path = paste("/Users/alaverre/Documents/HelmetedCurassowGenome/results/protein_acceleration", sp, "/", sep="/")

#######################################################################################
# Prepare Trees 
PreparedTree = paste0(path, sp, "_RERconverge_Trees.rds")

if (!file.exists(PreparedTree)){
  TreeFile = paste0(path, sp, "_CombineTrees.txt")
  
  message(c("Input File:", TreeFile))
  message("Reading Tree File...")
  
  Trees=readTrees(TreeFile)
  
  print("Saving RERconverge residuals...")
  saveRDS(Trees, file=PreparedTree)
  
}else{
  message("Tree preparation for RER already done!")
  Trees = readRDS(PreparedTree)
}

######################################################################################
# Estimate RER and save output
REREstimation = paste0(path, sp, "_RERconverge_residuals.rds")

if (!file.exists(REREstimation)){
  print("Estimating relative evolutionary rate...")
  RER = getAllResiduals(Trees, transform = "sqrt", weighted = T, scale = T, plot=F)
  
  print("Saving RERconverge residuals...")
  saveRDS(RER, file=REREstimation)
  
}else{
  message("RER estimation for RER already done!")
  RER = readRDS(REREstimation)
}
  
######################################################################################
### Defining phenotypes
Helmeted_birds = c("Casuarius_casuarius", "Bucorvus_abyssinicus", "Buceros_rhinoceros",
                   "Pauxi_pauxi", "Numida_meleagris",  "Anseranas_semipalmata", 
                   "Balearica_regulorum", "Anser_cygnoides")
Helmeted_squamates = c("Basiliscus_vittatus", "Chamaeleo_calyptratus")
Helmeted_all_species = c(Helmeted_birds, Helmeted_squamates)

Helmeted = get(paste0("Helmeted_", sp))
message("Considered Helmeted species:")
print(Helmeted)

# Generating paths for incomplete trees
pheno <- foreground2Paths(Helmeted, Trees, clade=clade) #terminal, ancestral or all

write(pheno, file=paste0(path, sp, "_RER_PhenoPath_", clade, ".txt"))

#######################################################################################
### Correlation
min.sp = ifelse(sp=="squamates", 5, 10)
min.pos = ifelse(sp=="squamates", 2, 5)
cor_all=correlateWithBinaryPhenotype(RER, pheno, min.sp=min.sp, min.pos=min.pos, weighted="auto")

ordered_cor = cor_all[order(cor_all$p.adj, cor_all$P),]

write.table(ordered_cor, file=paste0(path, sp, "_RER_correlations_clade_", clade, ".txt"), row.names=T, col.names=T, sep="\t",quote=F)

#######################################################################################

print("done!")

######################################################################################
