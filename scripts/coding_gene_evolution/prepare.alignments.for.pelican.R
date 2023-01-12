###########################################################################
###########################################################################

library(stringr)

pathAnnot="../../data/ensembl_annotations/Chicken/"

###########################################################################

for(spset in c("all_species", "without_chameleons")){

    if(spset=="all_species"){
        pathOrthoGroups="../../results/gene_families/OrthoFinder/all_species/iqtree/Results_Jan05/Phylogenetic_Hierarchical_Orthogroups/"
    }

    if(spset=="without_chameleons"){
        pathOrthoGroups="../../results/gene_families/OrthoFinder/without_chameleons/iqtree/Results_Sep21/Phylogenetic_Hierarchical_Orthogroups/"
    }

    for(dataset in c("birds", "all_species")){

        pathResults=paste("../../results/coding_gene_evolution/",spset, "/", dataset,sep="")

        pathCDS=paste(pathResults, "/CDS/",sep="")
        pathOutput=paste(pathResults, "/data_for_pelican/",sep="")

        ###########################################################################

        if(dir.exists(pathOutput)){
            print("path output already there")
        } else{
            system(paste("mkdr ",pathOutput))
        }

        ###########################################################################

        names=read.table(paste(pathAnnot, "GeneNames_Ensembl103.txt", sep=""), h=T, stringsAsFactors=F)
        rownames(names)=names[,1]

        ###########################################################################

        if(dataset == "all_species"){
            orthogroups=read.table(paste(pathOrthoGroups, "N0.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")
        }

        if(dataset == "birds"){
            orthogroups=read.table(paste(pathOrthoGroups, "N2.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")
        }

        ###########################################################################

        for(i in 1:nrow(orthogroups)){
            hog=orthogroups$HOG[i]
            og=orthogroups$OG[i]

            pathIn=paste(pathCDS, hog,"_",og,".aln.best.fas",sep="")

            if(file.exists(pathIn)){
                id.chicken=orthogroups[i,"Gallus_gallus"]
                id.chicken=unlist(strsplit(".", id.chicken))[1]

                if(id.chicken%in%rownames(names)){
                    name.chicken=names[id.chicken, 2]

                    if(length(grep("/",name.chicken))!=0){
                        name.chicken=str_replace(name.chicken, '/', '_')
                    }

                    if(!is.na(name.chicken) & name.chicken!=""){
                        system(paste("cp ",pathIn, " ",pathOutput, "/",id.chicken,"_",name.chicken,".fa",sep=""))
                    } else{
                        system(paste("cp ",pathIn, " ",pathOutput, "/",id.chicken,".fa",sep=""))
                    }
                } else{
                    system(paste("cp ",pathIn, " ",pathOutput, "/",hog,"_",og,".fa",sep=""))
                }
            }

            if(i%%1000==0){
                print(i)
            }
        }
    }
}

###########################################################################



