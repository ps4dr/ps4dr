library(data.table)
library(tidyr)
library(RecordLinkage)

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
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Drop NA values from each lists
drug.cor = lapply(drug.cor, function(x) x[!is.na(x)])
drug.cor = lapply(drug.cor, function(x) as.data.frame(x))
x= drug.cor[["Crohn's disease"]]

d = density(drug.cor[["Crohn's disease"]])

cor_score= x %>% drop_na()

cor_score=as.data.frame(do.call(rbind, drug.cor))
library(tidyr)
cor_score= x %>% drop_na()
#cor_score$drug = rownames(cor_score)
plot(cor_score$V1)

d = density(cor_score$Correlation.Score)
pdf(file = "./data/spia/anticor_CD.pdf")
plot(d,main = "Drug anti-correlation Scores for Crohn's Disease")
abline(v=median(cor_score$V1),col="blue")
text(median(cor_score$V1), 1, cex=1.1, round(median(cor_score$V1), digits=4))
dev.off()

write.csv(cor_score,file = "./data/spia/drug.cor10.cd.txt",quote = F)



# drug.dis.path =  merge(drug.path,dis.path,by="Name")
# colnames(drug.dis.path) = c("Pathways","Drug.Influence","Disease.Influence")
# drug.dis.path$anti.Correlation = ifelse(drug.dis.path$Drug.Influence == drug.dis.path$Disease.Influence,0,1)
# drug.dis.path$Drug.Influence = ifelse(drug.dis.path$Drug.Influence == "Activated",1,-1)
# drug.dis.path$Disease.Influence = ifelse(drug.dis.path$Disease.Influence == "Activated",1,-1)
# #w = pdist(drug.dis.path,indices.A = 1:nrow(drug.dis.path), indices.B = 1:nrow(drug.dis.path))
# drug.cor = cor(drug.dis.path$Drug.Influence, drug.dis.path$Disease.Influence)

load("./data/spia/drug.cor100.ad.RData")
names(drug.cor100)[i]
colnames(drug.cor100)[i]

gsub("\\..*","",names(drug.cor100)[i])



x=as.data.frame(do.call(rbind, drug.cor))
names(x) = "anticor"
x$anticor = as.numeric(x$anticor)
x= na.omit(x)
mean(x$anticor)
median(x$anticor)



