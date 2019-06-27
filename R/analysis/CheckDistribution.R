#' 4th script
#' summary:
#' This script calculates correation and dissimilarity scores for all drug-disease pathway pairs
#' plot distribution for all drug correlation scores for each disease 


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(RecordLinkage)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(Hmisc)))

#####################################################################
#TODO: Change to the directory where you cloned this repository
#~~~~~~~Using relative path~~~~~~~#
ensureFolder = function(folder) {
  if (! file.exists(folder)) {
    dir.create(folder)
  }
}

args = commandArgs(trailingOnly = TRUE)
resultsFolder = normalizePath(args[1])
ensureFolder(resultsFolder)
sprintf("Using results folder at %s", resultsFolder)

dataFolder = file.path(resultsFolder)
#####################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~load KEGG Drug SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load(file.path(dataFolder,"spia_output/spia_kegg_30Diseases_nopar.RData"))
# load(file.path(dataFolder,"spia_output/spia_kegg_47Diseases_drugGWAS_nopar.RData"))
load(file.path(dataFolder,"spia_output/spia_kegg_47Diseases_drugPdisease_nopar.RData"))
spia_drug_kegg = spia_kegg_47D
rm(spia_kegg_47D)
spia_drug_kegg = Filter(function(x) ! is.null(x), spia_drug_kegg) #delete empty df from list

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#
for (i in seq_along(spia_drug_kegg)) {
    spia_drug_kegg[[i]] = lapply(spia_drug_kegg[[i]], function(x) x[x$pNDE <= 0.05,])
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# remove empty drug frames
for (i in 1 : length(spia_drug_kegg)) {
    spia_drug_kegg[[i]] = spia_drug_kegg[[i]][lapply(spia_drug_kegg[[i]], length) > 1]
}


drug.path = vector('list', length(spia_drug_kegg)) # create list of lists
names(drug.path) = names(spia_drug_kegg)
for (i in seq_along(spia_drug_kegg)) {
    for (j in seq_along(spia_drug_kegg[[i]])) {
        drug.path[[i]][[j]] = spia_drug_kegg[[i]][[j]][, c(1, 11)] # use 2 for ID
        names(drug.path[[i]])[[j]] = names(spia_drug_kegg[[i]])[[j]]
    }
}

#' filter out diseases which has less than 1 drugs
for (i in seq_along(drug.path)) {
    drug.path[[i]] = Filter(function(x) ! dim(x)[1] == 0, drug.path[[i]])
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~load KEGG Disease SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load(file.path(dataFolder,"spia_output/spia_kegg_disease42.genes50_results.RData"))
load(file.path(dataFolder,"spia_output/spia_kegg_disease47.genes50_results.RData"))
# load(file.path(dataFolder,"spia_output/spia_kegg_degs_disease61.genes50_results.RData"))
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
    for (j in seq_along(spia_kegg)) {
        if (names(drug.path)[[i]] == names(spia_kegg)[j]) {
            dis.path[[i]] = spia_kegg[[j]][, c(1, 11)]
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
        drug.dis.path[[i]][[j]] = merge(dis.path[[i]], drug.path[[i]][[j]], by = "Name")
        names(drug.dis.path[[i]])[[j]] = names(drug.path[[i]])[[j]]
        names(drug.dis.path[[i]][[j]]) = c("Pathways", "Disease.Influence", "Drug.Influence")
        drug.dis.path[[i]][[j]]$Disease.Influence = ifelse(drug.dis.path[[i]][[j]]$Disease.Influence == "Activated", 1, - 1)
        drug.dis.path[[i]][[j]]$Drug.Influence = ifelse(drug.dis.path[[i]][[j]]$Drug.Influence == "Activated", 1, - 1)

        drug.Correlation[[i]][[j]] = round(cor(drug.dis.path[[i]][[j]]$Drug.Influence, drug.dis.path[[i]][[j]]$Disease.Influence), 2)
        names(drug.Correlation[[i]])[[j]] = names(drug.dis.path[[i]])[[j]]
        drug.Correlation[[i]][[j]] = as.data.frame(drug.Correlation[[i]][[j]])
        names(drug.Correlation[[i]][[j]]) = "Correlation.Score"
        drug.Correlation[[i]][[j]]$"Dissimilarity.Score" = round((sum(levenshteinDist(as.character(drug.dis.path[[i]][[j]]$Drug.Influence), as.character(drug.dis.path[[i]][[j]]$Disease.Influence))) * 100) / length(drug.dis.path[[i]][[j]]$Drug.Influence), 2)
        drug.Correlation[[i]][[j]]$"DrugPathway" = length(drug.dis.path[[i]][[j]]$Drug.Influence)
        drug.Correlation[[i]][[j]]$"DiseasePathway" = length(dis.path[[i]]$Name)
        drug.Correlation[[i]][[j]]$"affectedPathway" = round((drug.Correlation[[i]][[j]]$"DrugPathway" / drug.Correlation[[i]][[j]]$"DiseasePathway") * 100, 2)
        drug.Correlation[[i]][[j]]$"Disease" = names(dis.path)[[i]]
    }
}

##~~~~~create a single data.frame from all drugs in a disease~~~##

for (i in seq_along(drug.Correlation)) {
    drug.Correlation[[i]] = do.call(rbind, drug.Correlation[[i]])
    drug.Correlation[[i]]$Drug = rownames(drug.Correlation[[i]])
    drug.Correlation[[i]] = drug.Correlation[[i]][, c(7, 6, 1, 2, 3, 4, 5)]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug.Correlation = do.call(rbind, drug.Correlation) # list to df
drug.Correlation = filter(drug.Correlation, ! is.na(Correlation.Score)) # remove all rows with correlation score NA
drug.Correlation = split(drug.Correlation, drug.Correlation$Disease) # df to list
drug.Correlation = lapply(drug.Correlation, data.table) # make all df to data table

#~~~~Remove any disease with only positive correlation scores for all drugs ~~~#
for (i in seq_along(drug.Correlation)) {
    drug.Correlation[[i]] = drug.Correlation[[i]][! all(drug.Correlation[[i]]$Correlation.Score >= 0)]
}

drug.Correlation = Filter(function(x) dim(x)[1] >= 1, drug.Correlation) # remove empty disease lists


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' filter out drugs from each disease with correlationscore greater than -0.1
drug.shortlist = drug.Correlation
drug.shortlist = lapply(drug.shortlist, data.table)

for (i in seq_along(drug.shortlist)) {
    drug.shortlist[[i]] = drug.shortlist[[i]][drug.shortlist[[i]]$Correlation.Score <= - 0.4 & drug.shortlist[[i]]$'affectedPathway' >= 50]
}


#' filter out diseases which has less than 5 drugs
drug.shortlist = Filter(function(x) dim(x)[1] >= 1, drug.shortlist)
drug.shortlist.df = do.call(rbind, drug.shortlist)
drug.shortlist.df$Drug = tolower(drug.shortlist.df$Drug)
drug.shortlist.df$Drug = capitalize(drug.shortlist.df$Drug)
fwrite(drug.shortlist.df, file = file.path(dataFolder,"drug.shortlist.csv"))
save(drug.path,dis.path,drug.dis.path,drug.Correlation,drug.shortlist,file=file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))
load(file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))

drug.Correlation.df = do.call(rbind, drug.Correlation)
x = drug.Correlation.df[drug.Correlation.df$Correlation.Score >= 0.5 & drug.Correlation.df$affectedPathway >= 80]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Scatter Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugCor = do.call(rbind, drug.Correlation)
# drugCor$Drug = as.factor(drugCor$Drug)
drugCor = data.table(drugCor %>% drop_na())
drugCor = drugCor[order(drugCor$Correlation.Score, - drugCor$affectedPathway),]
# drugCor$Drug = as.character(drugCor$Drug)
drugCor$Drug <- tolower(drugCor$Drug)
drugCor$Drug = capitalize(drugCor$Drug)
# drugCor$Disease = capitalize(drugCor$Disease)

#' remove few extreme correlation scores

jpeg(file = file.path(dataFolder,"ScatterPlots_DrugPDisease_CorrelationScore.jpeg"), width = 2800, height = 1980, res = 200)
ggplot(drugCor, aes(x = affectedPathway, y = Correlation.Score, col = Disease)) +
    geom_point(size = 2, shape = 1) +
    labs(title = "Scatter Plots of Correlation Scores and affected pathways(%)") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Correlation Sore") +
    xlab("affected Pathway (%)")
dev.off()

#' additional
# jpeg(file=file.path(dataFolder,"ScatterPlots_highlighted_DrugPDisease_CorrelationScore.jpeg"), width=2800, height=1980, res=200)
# ggplot(drugCor,aes(x= affectedPathway,y=Correlation.Score, col=Disease)) + geom_point(size=3, shape=1, colour ="grey") + 
#   geom_point(data = cmap2, aes(x= DrugPathway,y=Correlation.Score, col=Disease), size=3) +
#   labs(title="Scatter Plots of Correlation Scores and affected pathways(%)") + 
#   theme(legend.position = "bottom",legend.title=element_text(size=10)) + 
#   theme(plot.title = element_text(hjust = 0.5)) + ylab("Correlation Sore") + xlab("affected Pathway (%)") +
#   geom_hline(aes(yintercept = 0)) +
#   geom_vline(aes(xintercept = 50))
# dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~~~~~Q-Q Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

jpeg(file = file.path(dataFolder,"qqPlots_DrugPDisease_CorrelationScore.jpeg"), width = 2800, height = 1980, res = 200)
p = list()
for (i in seq_along(drug.Correlation)) {
    p[[i]] = ggplot(drug.Correlation[[i]], aes(sample = Correlation.Score, col = Disease)) +
        stat_qq() +
        stat_qq_line() +
        labs(title = drug.Correlation[[i]]$Disease) +
        theme(legend.position = "none") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
}

grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]], p[[19]], p[[20]], p[[21]], p[[22]], p[[23]], p[[24]], p[[25]], p[[26]], p[[27]], p[[28]], p[[29]], p[[30]], ncol = 5, nrow = 6)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Density Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))

for (i in seq_along(drug.Correlation)) {
    for (j in seq_along(drug.Correlation[[i]])) {
        drug.Correlation[[i]]$Dissimilarity.Score = NULL
        drug.Correlation[[i]]$DrugPathway = NULL
        drug.Correlation[[i]]$DiseasePathway = NULL
        drug.Correlation[[i]]$Disease = NULL
        drug.Correlation[[i]]$affectedPathway = NULL
    }
    names(drug.Correlation[[i]])[[2]] = names(drug.Correlation)[[i]]
}

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
# drug.Correlation$`celiac disease` = NULL
density.score = lapply(drug.Correlation, melt)
density.score = do.call(rbind, density.score)
names(density.score) = c("Drug", "Diseases", "Correlation_Score")

# jpeg(file=file.path(dataFolder,"densityPlots_allDiseases_DrugGWAS.jpeg"), width=2800, height=1980, res=200)
# jpeg(file=file.path(dataFolder,"densityPlots_allDiseases_DrugPdisease_pval_30D.jpeg"), width=2800, height=1980, res=200)
ggplot(density.score, aes(x = Correlation_Score, fill = Diseases)) +
    geom_density(alpha = 0.25) +
    labs(title = "Distribution of Correlation Scores of all Drugs for each Diseases") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Density") +
    xlab("Correlation Score")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~Standization~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
load(file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))
#' Standization with scale function (Z-score normalization)
#' 
# drug.Correlation = do.call(rbind,drug.Correlation) # list to df
drug.Correlation = lapply(drug.Correlation, function(x) na.omit(x))
drug.Correlation = lapply(drug.Correlation, data.frame) # need to convert back to data frame since following operations do not work on data table

# scaling with z score transformation
drug.Cor.Scaled = drug.Correlation
for (i in seq_along(drug.Correlation)) {
    drug.Cor.Scaled[[i]] = lapply(drug.Correlation[[i]][3], function(x) scale(x, scale = TRUE))
    drug.Cor.Scaled[[i]]$durg = drug.Correlation[[i]]$Drug
    drug.Cor.Scaled[[i]]$affectedPathway = drug.Correlation[[i]]$affectedPathway
    #drug.Cor.Scaled[[i]]$disease = names(drug.Cor.Scaled)[[i]]
    drug.Cor.Scaled[[i]] = data.table(do.call(cbind, drug.Cor.Scaled[[i]]))
    drug.Cor.Scaled[[i]] = drug.Cor.Scaled[[i]][, c(2, 1, 3)]
    drug.Cor.Scaled[[i]]$V1 = as.numeric(drug.Cor.Scaled[[i]]$V1)
    names(drug.Cor.Scaled[[i]])[2] = names(drug.Correlation[[i]])[2]
}

drug.Cor.Scaled$`celiac disease` = NULL #it gives spike on the graph due to less variabllity in the data
drug.Cor.Scaled$`non-small cell lung carcinoma` = NULL

drug.Cor.Scaled = lapply(drug.Cor.Scaled, melt)
for (i in seq_along(drug.Cor.Scaled)) {
    drug.Cor.Scaled[[i]]$variable = names(drug.Cor.Scaled)[[i]]
}
drug.Cor.Scaled = do.call(rbind, drug.Cor.Scaled)
drug.Cor.Scaled = data.table(drug.Cor.Scaled %>% drop_na())
drug.Cor.Scaled = drug.Cor.Scaled[, c(1, 3, 4, 2)]
names(drug.Cor.Scaled) = c("Drug", "Disease", "Correlation.Score", "affectedPathway")
drug.Cor.Scaled$Disease = as.character(drug.Cor.Scaled$Disease)
drug.Cor.Scaled$affectedPathway = as.numeric(drug.Cor.Scaled$affectedPathway)
drug.Cor.Scaled = drug.Cor.Scaled[order(drug.Cor.Scaled$Correlation.Score, - drug.Cor.Scaled$affectedPathway),]
drug.Cor.Scaled$Drug <- tolower(drug.Cor.Scaled$Drug)
drug.Cor.Scaled$Drug = capitalize(drug.Cor.Scaled$Drug)
#drug.Cor.Scaled$Disease = capitalize(drug.Cor.Scaled$Disease)


# jpeg(file=file.path(dataFolder,"densityPlots_DrugPdisease_CorrelationScore_scaled.jpeg"), width=2800, height=1980, res=200)
ggplot(drug.Cor.Scaled, aes(x = Correlation.Score, fill = Disease)) +
    geom_density(alpha = 0.25) +
    labs(title = "Distribution of Correlation Scores of all Drugs for each Diseases") +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 12)) +
    ylab("Density") +
    xlab("Correlation.Score")
dev.off()


# jpeg(file=file.path(dataFolder,"ScatterPlots_DrugPDisease_scaled_CorrelationScore.jpeg"), width=2800, height=1980, res=200)
ggplot(drug.Cor.Scaled, aes(x = affectedPathway, y = Correlation.Score, col = Disease)) +
    geom_point(size = 2, shape = 1) +
    labs(title = "Scatter Plots of Correlation Scores and affected pathways(%)") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Correlation Sore") +
    xlab("affected Pathway (%)")
dev.off()
