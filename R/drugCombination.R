#' create pathway activities with the combination of two drugs
#' then use 2drugs combination pathways to calculate anti-correlation score 

comb1=list()
comb2=list()
set.seed(123)
for(i in seq_along(drug.path)){
  comb1[[i]] = combn(drug.path[[i]],2)
  comb2[[i]] = combn(names(drug.path[[i]]),2)
}


drug.comb.path = vector('list', length(drug.path)) # create list of lists
names(drug.comb.path) <- names(drug.path)

for (i in 1:(length(comb1))){
  for (j in 1:(length(comb1[[i]])/2)) {
    
    drug.comb.path[[i]][[j]] = data.table(merge(comb1[[i]][,j][1],comb1[[i]][,j][2],by="Name",all=TRUE))
    names(drug.comb.path[[i]])[[j]] = paste(comb2[[i]][,j][1],comb2[[i]][,j][2],sep = "_")
    # create a dummy column by copmaring disease and drug affects columns
    drug.comb.path[[i]][[j]][, Status.update:= ifelse(is.na(drug.comb.path[[i]][[j]]$Status.x),drug.comb.path[[i]][[j]]$Status.y,
                                                      ifelse(is.na(drug.comb.path[[i]][[j]]$Status.y),drug.comb.path[[i]][[j]]$Status.x,
                                                             ifelse(drug.comb.path[[i]][[j]]$Status.x!=drug.comb.path[[i]][[j]]$Status.y,"Neutral",drug.comb.path[[i]][[j]]$Status.y)))]
  }
}




drug.dis.path = vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.dis.path) <- names(spia_drug_kegg)
drug.cor = vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.cor) <- names(spia_drug_kegg)

for (i in seq_along(drug.path)) {
  for (j in seq_along(drug.path[[i]])) {
    drug.dis.path[[i]][[j]] = merge(dis.path[[i]],drug.path[[i]][[j]],by="Name")
    names(drug.dis.path[[i]])[[j]] <- names(spia_drug_kegg[[i]])[[j]]
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











