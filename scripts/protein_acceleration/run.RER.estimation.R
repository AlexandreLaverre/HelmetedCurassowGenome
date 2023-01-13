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

sp = "Duck" #args[1]
prefix = "scATAC" #args[2]

path = paste("/home/laverre/IPLOSS/results/NoncodingElements_Evolution", sp, prefix, "RERconverge/", sep="/") #"/ifb/data/mydatalocal/IPLOSS/"

#if (filter){prefix = paste0(prefix, "_nofilter")}

#######################################################################################
# Read trees for RERconverge and save output
TreeFile = paste0(path, prefix, "_CombineTrees.txt")

message(c("Input File:", TreeFile))
message("Reading Tree File...")
Trees=readTrees(TreeFile)

print("Saving RERconverge residuals...")
saveRDS(Trees, file=paste0(path, prefix, "_RERconverge_Trees.rds"))

######################################################################################
# Estimate RER and save output
print("Estimating relative evolutionary rate")
RER = getAllResiduals(Trees, transform = "sqrt", weighted = T, scale = T, plot=F)

print("Saving RERconverge residuals...")
saveRDS(RER, file=paste0(path, prefix, "_RERconverge_residuals.rds"))

print("done!")

######################################################################################
