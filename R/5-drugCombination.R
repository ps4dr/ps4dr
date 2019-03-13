#' 5th script
#' summary:
#' create pathway activities with the combination of two drugs
#' then use 2drugs combination pathways to calculate correlation score 


#####################################################################
#TODO: Change to the directory where you cloned this repository
setwd("/home/memon/projects/msdrp/")
#####################################################################

comb1=list()
comb2=list()
set.seed(123)
for(i in seq_along(drug.path)){
  comb1[[i]] = combn(drug.path[[i]],2)
  comb2[[i]] = combn(names(drug.path[[i]]),2)
}


drug.comb.path = vector('list', length(drug.path)) # create list of lists
names(drug.comb.path) = names(drug.path)

for (i in 1:(length(comb1))){
  for (j in 1:(length(comb1[[i]])/2)) {
    
    drug.comb.path[[i]][[j]] = data.table(merge(comb1[[i]][,j][1],comb1[[i]][,j][2],by="Name",all=TRUE))
    names(drug.comb.path[[i]])[[j]] = paste(comb2[[i]][,j][1],comb2[[i]][,j][2],sep = "_")
    # create a dummy column by copmaring disease and drug affects columns
    drug.comb.path[[i]][[j]][, Status.update:= 
                               ifelse(is.na(drug.comb.path[[i]][[j]]$Status.x),drug.comb.path[[i]][[j]]$Status.y,
                                      ifelse(is.na(drug.comb.path[[i]][[j]]$Status.y),drug.comb.path[[i]][[j]]$Status.x,
                                            ifelse(drug.comb.path[[i]][[j]]$Status.x!=drug.comb.path[[i]][[j]]$Status.y,"Neutral",drug.comb.path[[i]][[j]]$Status.y)))]
    drug.comb.path[[i]][[j]]$Status.x = NULL
    drug.comb.path[[i]][[j]]$Status.y = NULL
  }
}

save(drug.comb.path,file="./data/drug.comb.path.RData")
rm(comb1,comb2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug.dis.path = vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.dis.path) = names(spia_drug_kegg)
drug.cor = vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.cor) = names(spia_drug_kegg)

for (i in seq_along(drug.comb.path)) {
  for (j in seq_along(drug.comb.path[[i]])) {
    drug.dis.path[[i]][[j]] = merge(dis.path[[i]],drug.comb.path[[i]][[j]],by="Name")
    names(drug.dis.path[[i]])[[j]] = names(drug.comb.path[[i]])[[j]]
    names(drug.dis.path[[i]][[j]]) = c("Pathways","Disease.Influence","Drug.Influence")
    drug.dis.path[[i]][[j]]$Disease.Influence = ifelse(drug.dis.path[[i]][[j]]$Disease.Influence == "Activated",1,-1)
    drug.dis.path[[i]][[j]]$Drug.Influence = ifelse(drug.dis.path[[i]][[j]]$Drug.Influence == "Activated",1,-1)
    
    
    drug.cor[[i]][[j]] = cor(drug.dis.path[[i]][[j]]$Drug.Influence, drug.dis.path[[i]][[j]]$Disease.Influence)
    names(drug.cor[[i]])[[j]] = names(drug.dis.path[[i]])[[j]]
    drug.cor[[i]][[j]] = as.data.frame(drug.cor[[i]][[j]])
    names(drug.cor[[i]][[j]]) = "Correlation.Score"
    drug.cor[[i]][[j]]$"Dissimilarity.Score" = (sum(levenshteinDist(as.character(drug.dis.path[[i]][[j]]$Drug.Influence),as.character(drug.dis.path[[i]][[j]]$Disease.Influence))) * 100) / length(drug.dis.path[[i]][[j]]$Drug.Influence)
    drug.cor[[i]][[j]]$"DrugPathway" = length(drug.dis.path[[i]][[j]]$Drug.Influence)
    drug.cor[[i]][[j]]$"DiseasePathway" = length(dis.path[[i]]$Name)
  }
}



##~~~~~create a single data.frame from all drugs in a disease~~~##

#x=do.call(rbind,drug.cor[["Alzheimer's disease"]])
for (i in seq_along(drug.cor)) {
  drug.cor[[i]] = do.call(rbind,drug.cor[[i]])
  drug.cor[[i]]$Drug = rownames(drug.cor[[i]])
  drug.cor[[i]] = drug.cor[[i]][,c(5,1,2,3,4)]
}

save(drug.dis.path,drug.cor, file = "./data/drugCombination.correlationScore.RData")
load("./data/drugCombination.correlationScore.RData")
rm(drug.dis.path)
drugPair.Correlation = drug.cor
drug.Correlation = drug.cor
save(drugPair.Correlation,drug.Correlation, file = "./data/completeDrug.correlationScore.RData")
save(drug.Correlation, file = "./data/Drug.correlationScore1.RData")
save(drugPair.Correlation, file = "./data/Drug.correlationScore2.RData")

load("./data/completeDrug.correlationScore.RData")
drug.cor = drugPair.Correlation

for (i in seq_along(drug.cor)) {
  for (j in seq_along(drug.cor[[i]])) {
    drug.cor[[i]]$Dissimilarity.Score = NULL
    drug.cor[[i]]$DrugPathway = NULL
    drug.cor[[i]]$DiseasePathway = NULL
    
  }
  names(drug.cor[[i]])[[2]] = names(drug.cor)[[i]]
}

#following two disease mess up the plot 
drug.cor$epilepsy = NULL
drug.cor$`HIV infection` = NULL
density.score = lapply(drug.cor, melt)
density.score = do.call(rbind,density.score)
names(density.score) = c("Drug","Diseases","Correlation_Score")
                          
density.score = density.score %>% drop_na()
density.score = data.table(density.score)

jpeg(file="./data/densityPlots_allDiseases_drugPairs.jpeg", width=2800, height=1980, res=200)
ggplot(density.score,aes(x=Correlation_Score, fill=Diseases)) + geom_density(alpha=0.25) + 
  labs(title="Distribution of Correlation Scores of all Drug-Pairs for each Diseases")+ theme(legend.position = "bottom",legend.title=element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Density") + xlab("Correlation Score")
dev.off()



