#' 4th script
#' summary:
#' This script calculates correation and dissimilarity scores for all drug-disease pathway pairs
#' plot distribution for all drug correlation scores for each disease 


library(data.table)
library(tidyr)
library(dplyr)
library(RecordLinkage)
library(ggplot2)
library(gridExtra)
library(Hmisc)

#####################################################################
#TODO: Change to the directory where you cloned this repository
setwd("/home/memon/projects/msdrp/")
#####################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~load KEGG Drug SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load("./data/spia/spia_kegg_30Diseases_nopar.RData")
# load("./data/spia/spia_kegg_47Diseases_drugGWAS_nopar.RData")
load("./data/spia/spia_kegg_47Diseases_drugPdisease_nopar.RData")
spia_drug_kegg = spia_kegg_47D
rm(spia_kegg_47D)
spia_drug_kegg = Filter(function(x) !is.null(x), spia_drug_kegg) #delete empty df from list

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#
for (i in seq_along(spia_drug_kegg)){
  spia_drug_kegg[[i]] = lapply(spia_drug_kegg[[i]], function(x) x[x$pNDE <= 0.05,])
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# remove empty drug frames
for (i in 1:length(spia_drug_kegg)) {
  spia_drug_kegg[[i]] = spia_drug_kegg[[i]][lapply(spia_drug_kegg[[i]], length) > 1]
}


drug.path = vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.path) = names(spia_drug_kegg)
for(i in seq_along(spia_drug_kegg)){
  for(j in seq_along(spia_drug_kegg[[i]])){
    drug.path[[i]][[j]] = spia_drug_kegg[[i]][[j]][,c(1,11)] # use 2 for ID
    names(drug.path[[i]])[[j]] = names(spia_drug_kegg[[i]])[[j]]
  }
}

#' filter out diseases which has less than 1 drugs
for (i in seq_along(drug.path)) {
  drug.path[[i]] = Filter(function(x) !dim(x)[1] == 0, drug.path[[i]])
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~load KEGG Disease SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load("./data/spia/spia_kegg_disease42.genes50_results.RData")
load("./data/spia/spia_kegg_disease47.genes50_results.RData")
# load("./data/spia/spia_kegg_degs_disease61.genes50_results.RData")
# spia_kegg = spia_kegg_degs
# dis.path = lapply(spia_kegg, function(x) x[,c(1,11)])

#~~~~Remove any disease pathway with p.value (pNDE) >= 0.05 ~~~#
spia_kegg = lapply(spia_kegg, function(x) x[x$pNDE <= 0.05,])

#~~~Create disiase paths with only pathway name and their activity status~~~#
#~~~~Remove any other diseases which are not in drug.path~~~#
#~~~~so, both drug.path & dis.path are equivalent~~~~~~~~~~~#

dis.path = vector('list', length(drug.path)) # create list of lists
names(dis.path) = names(drug.path)

for (i in seq_along(drug.path)) {
  for (j in seq_along(spia_kegg)){
    if (names(drug.path)[[i]] == names(spia_kegg)[j]) {
      dis.path[[i]] = spia_kegg[[j]][,c(1,11)]
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~Calculate Correlation-Score~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug.dis.path = vector('list', length(drug.path)) # create list of lists
names(drug.dis.path) = names(drug.path)
drug.Correlation = vector('list', length(drug.path)) # create list of lists
names(drug.Correlation) = names(drug.path)

for (i in seq_along(drug.path)) {
  for (j in seq_along(drug.path[[i]])) {
    drug.dis.path[[i]][[j]] = merge(dis.path[[i]],drug.path[[i]][[j]],by="Name")
    names(drug.dis.path[[i]])[[j]] = names(drug.path[[i]])[[j]]
    names(drug.dis.path[[i]][[j]]) = c("Pathways","Disease.Influence","Drug.Influence")
    drug.dis.path[[i]][[j]]$Disease.Influence = ifelse(drug.dis.path[[i]][[j]]$Disease.Influence == "Activated",1,-1)
    drug.dis.path[[i]][[j]]$Drug.Influence = ifelse(drug.dis.path[[i]][[j]]$Drug.Influence == "Activated",1,-1)

    drug.Correlation[[i]][[j]] = round(cor(drug.dis.path[[i]][[j]]$Drug.Influence, drug.dis.path[[i]][[j]]$Disease.Influence),2)
    names(drug.Correlation[[i]])[[j]] = names(drug.dis.path[[i]])[[j]]
    drug.Correlation[[i]][[j]] = as.data.frame(drug.Correlation[[i]][[j]])
    names(drug.Correlation[[i]][[j]]) = "Correlation.Score"
    drug.Correlation[[i]][[j]]$"Dissimilarity.Score" = round((sum(levenshteinDist(as.character(drug.dis.path[[i]][[j]]$Drug.Influence),as.character(drug.dis.path[[i]][[j]]$Disease.Influence))) * 100) / length(drug.dis.path[[i]][[j]]$Drug.Influence),2)
    drug.Correlation[[i]][[j]]$"DrugPathway" = length(drug.dis.path[[i]][[j]]$Drug.Influence)
    drug.Correlation[[i]][[j]]$"DiseasePathway" = length(dis.path[[i]]$Name)
    drug.Correlation[[i]][[j]]$"affectedPathway" = round((drug.Correlation[[i]][[j]]$"DrugPathway" / drug.Correlation[[i]][[j]]$"DiseasePathway")*100,2)
    drug.Correlation[[i]][[j]]$"Disease" = names(dis.path)[[i]]
  }
}

##~~~~~create a single data.frame from all drugs in a disease~~~##

for (i in seq_along(drug.Correlation)) {
  drug.Correlation[[i]] = do.call(rbind,drug.Correlation[[i]])
  drug.Correlation[[i]]$Drug = rownames(drug.Correlation[[i]])
  drug.Correlation[[i]] = drug.Correlation[[i]][,c(7,6,1,2,3,4,5)]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug.Correlation = do.call(rbind,drug.Correlation) # list to df
drug.Correlation = filter(drug.Correlation, !is.na(Correlation.Score)) # remove all rows with correlation score NA
drug.Correlation = split(drug.Correlation,drug.Correlation$Disease) # df to list
drug.Correlation = lapply(drug.Correlation, data.table) # make all df to data table

#~~~~Remove any disease with only positive correlation scores for all drugs ~~~#
for (i in seq_along(drug.Correlation)) {
  drug.Correlation[[i]] = drug.Correlation[[i]][!all(drug.Correlation[[i]]$Correlation.Score >= 0)]
}

drug.Correlation = Filter(function(x) dim(x)[1] >= 1, drug.Correlation) # remove empty disease lists


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' filter out drugs from each disease with correlationscore greater than -0.1
drug.shortlist = drug.Correlation
drug.shortlist = lapply(drug.shortlist, data.table)

for (i in seq_along(drug.shortlist)) {
  drug.shortlist[[i]] = drug.shortlist[[i]][drug.shortlist[[i]]$Correlation.Score <= -0.4 & drug.shortlist[[i]]$'affectedPathway' >= 50]
}


#' filter out diseases which has less than 5 drugs
drug.shortlist = Filter(function(x) dim(x)[1] >= 1, drug.shortlist)
drug.shortlist.df = do.call(rbind,drug.shortlist)
drug.shortlist.df$Drug = tolower(drug.shortlist.df$Drug)
drug.shortlist.df$Drug = capitalize(drug.shortlist.df$Drug)
# fwrite(drug.shortlist.df, file = "./data/drug.shortlist.csv")
# save(drug.path,dis.path,drug.dis.path,drug.Correlation,drug.shortlist,file="./data/drugCorraltion.drugGWAS.RData")
# save(drug.path,dis.path,drug.dis.path,drug.Correlation,drug.shortlist,file="./data/drugCorraltion.drug_onlyDEGS.RData")
# save(drug.path,dis.path,drug.dis.path,drug.Correlation,drug.shortlist,file="./data/drugCorraltion-1.drugPdisease.RData")
# load("./data/drugCorraltion.drugGWAS.RData")
load("./data/drugCorraltion-1.drugPdisease.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Scatter Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugCor = do.call(rbind,drug.Correlation)
# drugCor$Drug = as.factor(drugCor$Drug)
drugCor = data.table(drugCor %>% drop_na())
drugCor = drugCor[order(drugCor$Correlation.Score,-drugCor$affectedPathway),]
# drugCor$Drug = as.character(drugCor$Drug)
drugCor$Drug <- tolower(drugCor$Drug)
drugCor$Drug = capitalize(drugCor$Drug)
# drugCor$Disease = capitalize(drugCor$Disease)

#' remove few extreme correlation scores

jpeg(file="./data/ScatterPlots_DrugPDisease_CorrelationScore.jpeg", width=2800, height=1980, res=200)
ggplot(drugCor,aes(x= affectedPathway,y=Correlation.Score, col=Disease)) + geom_point(size=2, shape=1) +
  labs(title="Scatter Plots of Correlation Scores and affected pathways(%)") + 
  theme(legend.position = "bottom",legend.title=element_text(size=10)) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Correlation Sore") + xlab("affected Pathway (%)")
dev.off()

jpeg(file="./data/ScatterPlots_highlighted_DrugPDisease_CorrelationScore.jpeg", width=2800, height=1980, res=200)
ggplot(drugCor,aes(x= affectedPathway,y=Correlation.Score, col=Disease)) + geom_point(size=3, shape=1, colour ="grey") + 
  geom_point(data = cmap2, aes(x= DrugPathway,y=Correlation.Score, col=Disease), size=3) +
  labs(title="Scatter Plots of Correlation Scores and affected pathways(%)") + 
  theme(legend.position = "bottom",legend.title=element_text(size=10)) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Correlation Sore") + xlab("affected Pathway (%)") +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 50))
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~~~~~Q-Q Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

jpeg(file="./data/qqPlots_DrugPDisease_CorrelationScore.jpeg", width=2800, height=1980, res=200)
p = list()
for (i in seq_along(drug.Correlation)) {
  p[[i]] = ggplot(drug.Correlation[[i]],aes(sample=Correlation.Score, col=Disease)) + stat_qq() + stat_qq_line() +
    labs(title=drug.Correlation[[i]]$Disease) + theme(legend.position = "none")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank())
}

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]], ncol=5, nrow=6)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in seq_along(drug.Correlation)) {
  for (j in seq_along(drug.Correlation[[i]])) {
    drug.Correlation[[i]]$Dissimilarity.Score = NULL
    drug.Correlation[[i]]$DrugPathway = NULL
    drug.Correlation[[i]]$DiseasePathway = NULL
    drug.Correlation[[i]]$Disease = NULL
    # drug.Correlation[[i]]$affectedPathway = NULL
    
  }
  names(drug.Correlation[[i]])[[2]] = names(drug.Correlation)[[i]]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~Standization~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Standization with scale function (Z-score normalization)
#' 
# drug.Correlation = do.call(rbind,drug.Correlation) # list to df
drug.Correlation = lapply(drug.Correlation, function(x) na.omit(x))
drug.Correlation = lapply(drug.Correlation, data.frame) # need to convert back to data frame since following operations do not work on data table

# scaling with z score transformation
drug.Cor.Scaled = drug.Correlation
for (i in seq_along(drug.Correlation)) {
  drug.Cor.Scaled[[i]] = lapply(drug.Correlation[[i]][2],function(x) scale(x, scale=TRUE))
  drug.Cor.Scaled[[i]]$durg = drug.Correlation[[i]]$Drug
  drug.Cor.Scaled[[i]]$affectedPathway = drug.Correlation[[i]]$affectedPathway
  #drug.Cor.Scaled[[i]]$disease = names(drug.Cor.Scaled)[[i]]
  drug.Cor.Scaled[[i]] = data.table(do.call(cbind,drug.Cor.Scaled[[i]]))
  drug.Cor.Scaled[[i]] = drug.Cor.Scaled[[i]][,c(2,1,3)]
  drug.Cor.Scaled[[i]]$V1 = as.numeric(drug.Cor.Scaled[[i]]$V1)
  names(drug.Cor.Scaled[[i]])[2] = names(drug.Correlation[[i]])[2]
}

drug.Cor.Scaled$`celiac disease` = NULL #it gives spike on the graph due to less variabllity in the data
drug.Cor.Scaled$`non-small cell lung carcinoma` = NULL

drug.Cor.Scaled = lapply(drug.Cor.Scaled, melt)
for (i in seq_along(drug.Cor.Scaled)) {
  drug.Cor.Scaled[[i]]$variable = names(drug.Cor.Scaled)[[i]]
}
drug.Cor.Scaled = do.call(rbind,drug.Cor.Scaled)
drug.Cor.Scaled = data.table(drug.Cor.Scaled %>% drop_na())
drug.Cor.Scaled = drug.Cor.Scaled[,c(1,3,4,2)]
names(drug.Cor.Scaled) = c("Drug","Disease","Correlation.Score","affectedPathway")
drug.Cor.Scaled$Disease = as.character(drug.Cor.Scaled$Disease)
drug.Cor.Scaled$affectedPathway = as.numeric(drug.Cor.Scaled$affectedPathway)
drug.Cor.Scaled = drug.Cor.Scaled[order(drug.Cor.Scaled$Correlation.Score,-drug.Cor.Scaled$affectedPathway),]
drug.Cor.Scaled$Drug <- tolower(drug.Cor.Scaled$Drug)
drug.Cor.Scaled$Drug = capitalize(drug.Cor.Scaled$Drug)
#drug.Cor.Scaled$Disease = capitalize(drug.Cor.Scaled$Disease)


jpeg(file="./data/densityPlots_DrugPdisease_CorrelationScore_scaled.jpeg", width=2800, height=1980, res=200)
ggplot(drug.Cor.Scaled,aes(x=Correlation.Score, fill=Disease)) + geom_density(alpha=0.25) + 
  labs(title="Distribution of Correlation Scores of all Drugs for each Diseases")+ theme(legend.position = "bottom",legend.title=element_text(size=10)) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Density") + xlab("Correlation.Score")
dev.off()


jpeg(file="./data/ScatterPlots_DrugPDisease_scaled_CorrelationScore.jpeg", width=2800, height=1980, res=200)
ggplot(drug.Cor.Scaled,aes(x= affectedPathway,y=Correlation.Score, col=Disease)) + geom_point(size=2, shape=1) +
  labs(title="Scatter Plots of Correlation Scores and affected pathways(%)") + 
  theme(legend.position = "bottom",legend.title=element_text(size=10)) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Correlation Sore") + xlab("affected Pathway (%)")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~Drug mapping to Diseases~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugsA = unique(drug.Cor.Scaled$Drug)

#' this drug-indications table is downloaded from ChEMBL https://www.ebi.ac.uk/chembl/drug/indications on 15 March, 2019
chembl = fread("/home/memon/projects/msdrp/data/chembl_indications-19_10_36_41.txt")
chembl = chembl[,c(1,2,5:9)]
chembl = chembl[MAX_PHASE_FOR_IND ==1 |MAX_PHASE_FOR_IND ==2 | MAX_PHASE_FOR_IND ==3 | MAX_PHASE_FOR_IND == 4] 
#chembl = chembl[MAX_PHASE_FOR_IND ==3 | MAX_PHASE_FOR_IND == 4] 
chembl2 = chembl[,c(2,6,7)]
#chembl2$MOLECULE_NAME = toupper(chembl2$MOLECULE_NAME)
# drugsB = unique(toupper(chembl2$MOLECULE_NAME))
drugsB = unique(chembl2$MOLECULE_NAME)

DrugCommon = as.data.frame(intersect(drugsB,drugsA))
names(DrugCommon) = "Drug"
length(intersect(unique(chembl2$EFO_NAME),unique(drug.Cor.Scaled$Disease)))

#' update with drugs which are present in chemble table 
drug.Cor.Scaled = merge(drug.Cor.Scaled,DrugCommon,by="Drug") 

#' Clinically trialled drugs that appears in our analysis
cmap = merge(drug.Cor.Scaled, chembl2, by.x = c("Drug","Disease"),by.y=c("MOLECULE_NAME","EFO_NAME"))
# cmap2 = merge(drugCor, chembl2, by.x = c("Drug","Disease"),by.y=c("MOLECULE_NAME","EFO_NAME"))
names(cmap)[4] = "phase"
length(unique(cmap$Drug))
length(unique(cmap$Disease))

#' compute 10th percentile scores for each diseases
#' 
library(dplyr)
library(purrr)
qn10 = cmap %>% 
  dplyr::group_by(Disease) %>% 
  dplyr::summarize(quants = quantile(Correlation.Score, probs = 0.1))
qn10 = split(qn10, qn10$Disease)
cmap2 = split(cmap,cmap$Disease)

#' get all approved drugs that appear in 10th percentile
for (i in seq_along(cmap2)) {
  cmap2[[i]] = cmap2[[i]][Correlation.Score <= qn10[[i]][["quants"]]]
}
cmap3 = do.call(rbind,cmap2)

#' get any drugs that appear in 10th percentile
disease21 = as.data.frame(unique(cmap3$Disease))
names(disease21) = "Disease"
drug.Cor.Scaled2 = data.table(merge(disease21, drug.Cor.Scaled,by= "Disease"))

qn10all = drug.Cor.Scaled2 %>% 
  dplyr::group_by(Disease) %>% 
  dplyr::summarize(quants = quantile(Correlation.Score, probs = 0.1))
qn10all = split(qn10all, qn10all$Disease)

drug.Cor.Scaled2 = split(drug.Cor.Scaled2,drug.Cor.Scaled2$Disease)
for (i in seq_along(drug.Cor.Scaled2)) {
  drug.Cor.Scaled2[[i]] = drug.Cor.Scaled2[[i]][Correlation.Score <= qn10all[[i]][["quants"]]]
}
drug.Cor.Scaled2 = do.call(rbind,drug.Cor.Scaled2)


#~~~~~~~#
library(dplyr)
topDrugs = count(cmap3,Disease)
allDrugs = count(drug.Cor.Scaled,Disease)
Drugs10thParc = count(drug.Cor.Scaled2,Disease)
CombinedDrugs = merge(allDrugs,topDrugs,by="Disease")
CombinedDrugs = merge(CombinedDrugs,Drugs10thParc,by="Disease")
names(CombinedDrugs) = c("Disease","allDrugs","FDADrugs","Drugs10thP")
CombinedDrugs$percentage1 = round((CombinedDrugs$FDADrugs/CombinedDrugs$allDrugs)*100,2)
CombinedDrugs$percentage2 = round((CombinedDrugs$FDADrugs/CombinedDrugs$Drugs10thP)*100,2)
CombinedDrugs = CombinedDrugs[order(-CombinedDrugs$percentage2),]

write.csv(CombinedDrugs,"./data/CombinedDrugs.csv",quote=FALSE,row.names=FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Density Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




# following drug(s) creates anomaly in the graph due to smaller density distribution in DrugGWAS data
# drug.Correlation$`chronic obstructive pulmonary disease` = NULL
# 
# # following drug(s) creates anomaly in the graph due to smaller density distribution in DrugGWAS data
# # drug.Correlation$`HIV infection` = NULL
# # drug.Correlation$`celiac disease` = NULL
# # drug.Correlation$`mucocutaneous lymph node syndrome` = NULL
# drug.Correlation$`central nervous system cancer` = NULL
# drug.Correlation$`Pick disease` = NULL
# drug.Correlation$`Parkinson's disease` = NULL
# drug.Correlation$`mucocutaneous lymph node syndrome` = NULL
# drug.Correlation$`non-small cell lung carcinoma` = NULL
# drug.Correlation$`non-small cell lung carcinoma` = NULL

density.score = lapply(drug.Correlation, melt)
density.score = do.call(rbind,density.score)
names(density.score) = c("Drug","Diseases","Correlation_Score")

# jpeg(file="./data/densityPlots_allDiseases_DrugGWAS.jpeg", width=2800, height=1980, res=200)
# jpeg(file="./data/densityPlots_allDiseases_DrugPdisease_pval.jpeg", width=2800, height=1980, res=200)
ggplot(density.score,aes(x=Correlation_Score, fill=Diseases)) + geom_density(alpha=0.25) + 
  labs(title="Distribution of Correlation Scores of all Drugs for each Diseases")+ theme(legend.position = "bottom",legend.title=element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Density") + xlab("Correlation Score")
dev.off() 





##~~~extra

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~Drug mapping to Diseases~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
drugCor = do.call(rbind,drug.Correlation)
drugsA = unique(drugCor$Drug)

#' this drug disease relationships are fetched from OpenTargets
load("./data/drug2disease.RData")
d2d = drug2disease[!phase == 0]
d2d = unique(d2d[,c(4,2)])
drugsB = unique(d2d$chembl.name)

length(intersect(drugsA,drugsB))

#' this drug-indications table is downloaded from ChEMBL https://www.ebi.ac.uk/chembl/drug/indications on 15 March, 2019
chembl = fread("/home/memon/projects/msdrp/data/chembl_indications-19_10_36_41.txt")
chembl = chembl[,c(1,2,5:9)]
chembl = chembl[MAX_PHASE_FOR_IND ==1 |MAX_PHASE_FOR_IND ==2 | MAX_PHASE_FOR_IND ==3 | MAX_PHASE_FOR_IND == 4] 
chembl2 = chembl[,c(2,6,7)]
chembl2$MOLECULE_NAME = toupper(chembl2$MOLECULE_NAME)
drugsC = unique(toupper(chembl$MOLECULE_NAME))

length(intersect(drugsA,drugsC))

cmap = merge(drugCor, chembl2, by.x = c("Drug","Disease"),by.y=c("MOLECULE_NAME","EFO_NAME"))
#cmap2 = merge(drugCor, d2d, by.x = c("Drug","Disease"),by.y=c("chembl.name","efo.term"))
names(cmap)[8] = "phase"

#~~~~~~~#

length(unique(cmap$Drug))
length(unique(cmap$Disease))

dis.list = split(cmap, cmap$Disease)
dis.list = lapply(dis.list, data.table)

min.phase = 3
fTest = list()
for (i in seq_along(dis.list)) {
  a = nrow(dis.list[[i]][Correlation.Score < 0 & phase >= min.phase])
  b = nrow(dis.list[[i]][Correlation.Score > 0 & phase >= min.phase])
  c = nrow(dis.list[[i]][Correlation.Score < 0 & phase < min.phase])
  d = nrow(dis.list[[i]][Correlation.Score > 0 & phase < min.phase])
  fisher = fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  fTest[[i]] = data.table(odds.ratio = fisher$estimate, p.value = fisher$p.value,approved_negative = a,approved_positive = b)
  names(fTest)[[i]] = names(dis.list)[[i]]
  fTest[[i]]$Disease = names(dis.list)[[i]]
}

fisherall = do.call(rbind,fTest,envir = parent.frame())
fisherall = fisherall[,c(5,3,4,1,2)]


test=fread("/home/memon/projects/msdrp/data/chembl_drugtargets-19_14_39_48.txt")
test=fread("/home/memon/projects/msdrp/data/drugbank_dtis.tsv")
