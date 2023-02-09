###########################################################################
###########################################################################

library(stringr)

pathAnnot="../../data/ensembl_annotations/"

###########################################################################

for(spset in c("all_species")){ ## other option is "without_chameleons"

    if(spset=="all_species"){
        pathOrthoGroups="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/Phylogenetic_Hierarchical_Orthogroups/"
    }

    if(spset=="without_chameleons"){
        pathOrthoGroups="../../results/gene_families/OrthoFinder/without_chameleons/iqtree/Results_Sep21/Phylogenetic_Hierarchical_Orthogroups/"
    }

    for(dataset in c("birds", "all_species", "squamates")){

        pathResults=paste("../../results/coding_gene_evolution/",spset, "/", dataset,sep="")

        pathCDS=paste(pathResults, "/CDS/",sep="")
        pathOutput=paste(pathResults, "/data_for_pelican/",sep="")

        ###########################################################################

        if(dir.exists(pathOutput)){
            print("path output already there")
        } else{
            system(paste("mkdir ",pathOutput))
        }

        ###########################################################################

        if(dataset == "all_species"){
            orthogroups=read.table(paste(pathOrthoGroups, "N0.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")
        }

        if(dataset == "birds"){
            orthogroups=read.table(paste(pathOrthoGroups, "N2.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")
        }

        if(dataset == "squamates"){
            orthogroups=read.table(paste(pathOrthoGroups, "N1.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")
        }

        ###########################################################################

        for(i in 1:nrow(orthogroups)){
            hog=orthogroups$HOG[i]
            og=orthogroups$OG[i]

            pathIn=paste(pathCDS, hog,"_",og,".aln.best.fas",sep="")
            system(paste("cp ",pathIn, " ",pathOutput, "/",hog,"_",og,".fa",sep=""))

            if(i%%1000==0){
                print(i)
            }
        }
    }
}

###########################################################################



