library(data.table)
library(tidyr)
library(RecordLinkage)
library(ggplot2)

setwd("/home/memon/projects/msdrp/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~load KEGG Drug SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
load("./data/spia/spia_kegg_30Diseases_nopar.RData")
spia_drug_kegg = spia_kegg_30D
spia_drug_kegg = Filter(function(x) !is.null(x), spia_drug_kegg) #delete empty df from list

drug.path <- vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.path) <- names(spia_drug_kegg)
for(i in seq_along(spia_drug_kegg)){
  for(j in seq_along(spia_drug_kegg[[i]])){
    drug.path[[i]][[j]] <- spia_drug_kegg[[i]][[j]][,c(1,11)]
    names(drug.path[[i]])[[j]] <- names(spia_drug_kegg[[i]])[[j]]
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~load KEGG Disease SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("./data/spia/spia_kegg_disease42.genes50_results.RData")
# dis.path = lapply(spia_kegg, function(x) x[,c(1,11)])

#~~~~Remove any other diseases which are not in drug.path~~~#
#~~~~so, both drug.path & dis.path are equivalent~~~~~~~~~~~#

dis.path = vector('list', length(spia_drug_kegg)) # create list of lists
names(dis.path) <- names(spia_drug_kegg)

for (i in seq_along(spia_drug_kegg)) {
  for (j in seq_along(spia_kegg)){
    if (names(spia_drug_kegg)[[i]] == names(spia_kegg)[j]) {
      dis.path[[i]] = spia_kegg[[j]][,c(1,11)]
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~Calculate Correlation-Score~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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


##~~~~~create a single data.frame from all drugs in a disease~~~##

#x=do.call(rbind,drug.cor[["Alzheimer's disease"]])
for (i in seq_along(drug.cor)) {
  drug.cor[[i]] = do.call(rbind,drug.cor[[i]])
  drug.cor[[i]]$Drug = rownames(drug.cor[[i]])
  drug.cor[[i]] = drug.cor[[i]][,c(5,1,2,3,4)]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Density Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug.cor1 = drug.cor

for (i in seq_along(drug.cor)) {
  for (j in seq_along(drug.cor[[i]])) {
    drug.cor[[i]]$Dissimilarity.Score = NULL
    drug.cor[[i]]$DrugPathway = NULL
    drug.cor[[i]]$DiseasePathway = NULL
    
  }
  names(drug.cor[[i]])[[2]] = names(drug.cor)[[i]]
}

drug.cor$`HIV infection` = NULL
density.score = lapply(drug.cor, melt)
density.score = do.call(rbind,density.score)
names(density.score) <- c("Drug","Diseases","Correlation_Score")

jpeg(file="./data/densityPlots_allDiseases.jpeg", width=2800, height=1980, res=200)
ggplot(density.score,aes(x=Correlation_Score, fill=Diseases)) + geom_density(alpha=0.25) + 
  labs(title="Distribution of Correlation Scores of all Drugs for each Diseases")+ theme(legend.position = "bottom",legend.title=element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Density") + xlab("Correlation Score")
dev.off()


density.score = data.table(density.score)
density.score1 = density.score[variable == "hepatocellular carcinoma"]
density.score1= density.score1 %>% drop_na()

q5 <- quantile(density.score1$value,.05)
q95 <- quantile(density.score1$value,.95)


ggplot(density.score1,aes(x=value, fill=variable)) + geom_density(alpha=.25)




drug.cor = drug.cor1
cor_score= drug.cor %>% drop_na()
d = density(cor_score$Correlation.Score)
plot(d,main = "Drug anti-correlation Scores for Crohn's Disease")
abline(v=median(cor_score$V1),col="blue")
text(median(cor_score$V1), 1, cex=1.1, round(median(cor_score$V1), digits=4))


