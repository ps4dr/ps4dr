library(data.table)
#library(pdist)

setwd("/home/memon/projects/msdrp/")

#####_____load KEGG disease SPIA results_____#####
# load("~/projects/drugrep/drug2path/data/spia/spia_kegg_results.RData")
load("./data/spia/spia_kegg_disease42.genes50_results.RData")
# dis.path = data.table(spia_kegg[["EFO_0000400"]][,c(1,11)])
#dis.path = data.table(spia_kegg[["EFO_0000249"]][,c(1,11)])
# dis.path = data.table(spia_kegg[["EFO_0002508"]][,c(1,11)])
dis.path = data.table(spia_kegg[["Crohn's disease"]][,c(1,11)])


#####_____load KEGG Drug SPIA results_____#####
# load("~/projects/drugrep/drug2path/data/spia/spia_test_nopar_lfcdis1_try1.RData")
# spia_drug_kegg = spia_test$EFO_0000249
# spia_drug_kegg = spia_test[["EFO_0002508"]]
# spia_drug_kegg = spia_test[["EFO_0000400"]]
load("./data/spia/spia_kegg_30Diseases_nopar.RData")
spia_drug_kegg = spia_kegg_30D[["Crohn's disease"]]

spia_drug_kegg = Filter(function(x) !is.null(x), spia_drug_kegg) #delete empty df from list

drug.path <- list()
for (i in 1:length(spia_drug_kegg)) {
  drug.path[[i]] = spia_drug_kegg[[i]][,c(1,11)]
  names(drug.path)[i] = names(spia_drug_kegg)[i]
}

drug.dis.path = list()
drug.cor = list()
for (i in 1:length(drug.path)) {
  drug.dis.path[[i]] = merge(drug.path[[i]],dis.path,by="Name")
  names(drug.dis.path)[i] = names(spia_drug_kegg)[i]
  names(drug.dis.path[[i]]) = c("Pathways","Drug.Influence","Disease.Influence")
  drug.dis.path[[i]]$Drug.Influence = ifelse(drug.dis.path[[i]]$Drug.Influence == "Activated",1,-1)
  drug.dis.path[[i]]$Disease.Influence = ifelse(drug.dis.path[[i]]$Disease.Influence == "Activated",1,-1)
  drug.cor[[i]] = cor(drug.dis.path[[i]]$Drug.Influence, drug.dis.path[[i]]$Disease.Influence)
  names(drug.cor)[i] = gsub("\\..*","",names(drug.dis.path)[i])
}

cor_score=as.data.frame(do.call(rbind, drug.cor))
library(tidyr)
cor_score= cor_score %>% drop_na()
#cor_score$drug = rownames(cor_score)
plot(cor_score$V1)

d = density(cor_score$V1)
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


