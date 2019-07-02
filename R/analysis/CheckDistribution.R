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
suppressWarnings(suppressMessages(library(purrr)))

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

load(file.path(dataFolder,"spia_output/spia_kegg_drug_Diseases_drugPdisease_nopar.RData"))
spia_kegg_drug = Filter(function(x) ! is.null(x), spia_kegg_drug) #delete empty df from list

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#
for (i in seq_along(spia_kegg_drug)) {
    spia_kegg_drug[[i]] = lapply(spia_kegg_drug[[i]], function(x) x[x$pNDE <= 0.05,])
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# remove empty drug frames
for (i in 1 : length(spia_kegg_drug)) {
    spia_kegg_drug[[i]] = spia_kegg_drug[[i]][lapply(spia_kegg_drug[[i]], length) > 1]
}


drug_path_temp = vector('list', length(spia_kegg_drug)) # create list of lists
names(drug_path_temp) = names(spia_kegg_drug)
for (i in seq_along(spia_kegg_drug)) {
    for (j in seq_along(spia_kegg_drug[[i]])) {
        drug_path_temp[[i]][[j]] = spia_kegg_drug[[i]][[j]][, c(1, 11)] # use 2 for ID
        names(drug_path_temp[[i]])[[j]] = names(spia_kegg_drug[[i]])[[j]]
    }
}

#' filter out diseases which has less than 1 drugs
for (i in seq_along(drug_path_temp)) {
    drug_path_temp[[i]] = Filter(function(x) ! dim(x)[1] == 0, drug_path_temp[[i]])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~load KEGG Disease SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_output/spia_kegg_diseaseGenes50.RData"))

#~~~~Remove any disease pathway with p.value (pNDE) >= 0.05 ~~~#
spia_kegg = lapply(spia_kegg, function(x) x[x$pNDE <= 0.05,])
spia_kegg = discard(spia_kegg, ~ all(is.na(.x))) # remove empty diseases


#~~~Create disiase paths with only pathway name and their activity status~~~#
#~~~~Remove any other diseases which are not in drug_path_temp~~~#
#~~~~so, both drug_path_temp & dis_path are equivalent~~~~~~~~~~~#

dis_path = vector('list', length(drug_path_temp)) # create list of lists
names(dis_path) = names(drug_path_temp)

for (i in seq_along(drug_path_temp)) {
    for (j in seq_along(spia_kegg)) {
        if (names(drug_path_temp)[[i]] == names(spia_kegg)[j]) {
            dis_path[[i]] = spia_kegg[[j]][, c(1, 11)]
        }
    }
}

# remove empty disease frames
dis_path = discard(dis_path, ~ all(is.na(.x))) # remove empty diseases


#' Now, remove those diseases from drug_path_temp which are not in dis_path

drug_path = vector('list', length(dis_path)) # create list of lists
names(drug_path) = names(dis_path)

for (i in 1 : length(drug_path_temp)) {
  if (names(drug_path_temp)[[i]] %in% names(spia_kegg)) {
    drug_path[[i]] = drug_path_temp[[i]]
    names(drug_path)[[i]] = names(drug_path_temp)[[i]]
  }
}

drug_path = discard(drug_path, ~ all(is.na(.x))) # remove empty diseases

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~Calculate Correlation-Score~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_dis_path = vector('list', length(drug_path)) # create list of lists
names(drug_dis_path) = names(drug_path)
drug_correlation = vector('list', length(drug_path)) # create list of lists
names(drug_correlation) = names(drug_path)

for (i in seq_along(drug_path)) {
    for (j in seq_along(drug_path[[i]])) {
        drug_dis_path[[i]][[j]] = merge(dis_path[[i]], drug_path[[i]][[j]], by = "Name")
        names(drug_dis_path[[i]])[[j]] = names(drug_path[[i]])[[j]]
        names(drug_dis_path[[i]][[j]]) = c("Pathways", "Disease.Influence", "Drug.Influence")
        drug_dis_path[[i]][[j]]$Disease.Influence = ifelse(drug_dis_path[[i]][[j]]$Disease.Influence == "Activated", 1, - 1)
        drug_dis_path[[i]][[j]]$Drug.Influence = ifelse(drug_dis_path[[i]][[j]]$Drug.Influence == "Activated", 1, - 1)

        drug_correlation[[i]][[j]] = round(cor(drug_dis_path[[i]][[j]]$Drug.Influence, drug_dis_path[[i]][[j]]$Disease.Influence), 2)
        names(drug_correlation[[i]])[[j]] = names(drug_dis_path[[i]])[[j]]
        drug_correlation[[i]][[j]] = as.data.frame(drug_correlation[[i]][[j]])
        names(drug_correlation[[i]][[j]]) = "Correlation.Score"
        drug_correlation[[i]][[j]]$"Dissimilarity.Score" = round((sum(levenshteinDist(as.character(drug_dis_path[[i]][[j]]$Drug.Influence), as.character(drug_dis_path[[i]][[j]]$Disease.Influence))) * 100) / length(drug_dis_path[[i]][[j]]$Drug.Influence), 2)
        drug_correlation[[i]][[j]]$"DrugPathway" = length(drug_dis_path[[i]][[j]]$Drug.Influence)
        drug_correlation[[i]][[j]]$"DiseasePathway" = length(dis_path[[i]]$Name)
        drug_correlation[[i]][[j]]$"affectedPathway" = round((drug_correlation[[i]][[j]]$"DrugPathway" / drug_correlation[[i]][[j]]$"DiseasePathway") * 100, 2)
        drug_correlation[[i]][[j]]$"Disease" = names(dis_path)[[i]]
    }
}

##~~~~~create a single data.frame from all drugs in a disease~~~##

for (i in seq_along(drug_correlation)) {
    drug_correlation[[i]] = do.call(rbind, drug_correlation[[i]])
    drug_correlation[[i]]$Drug = rownames(drug_correlation[[i]])
    drug_correlation[[i]] = drug_correlation[[i]][, c(7, 6, 1, 2, 3, 4, 5)]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_correlation = do.call(rbind, drug_correlation) # list to df
drug_correlation = filter(drug_correlation, ! is.na(Correlation.Score)) # remove all rows with correlation score NA
drug_correlation = split(drug_correlation, drug_correlation$Disease) # df to list
drug_correlation = lapply(drug_correlation, data.table) # make all df to data table

#~~~~Remove any disease with only positive correlation scores for all drugs ~~~#
for (i in seq_along(drug_correlation)) {
    drug_correlation[[i]] = drug_correlation[[i]][! all(drug_correlation[[i]]$Correlation.Score >= 0)]
}

drug_correlation = Filter(function(x) dim(x)[1] >= 1, drug_correlation) # remove empty disease lists


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' filter out drugs from each disease with correlationscore greater than -0.1
drug_shortlist = drug_correlation
drug_shortlist = lapply(drug_shortlist, data.table)

for (i in seq_along(drug_shortlist)) {
    drug_shortlist[[i]] = drug_shortlist[[i]][drug_shortlist[[i]]$Correlation.Score <= - 0.4 & drug_shortlist[[i]]$'affectedPathway' >= 50]
}


#' filter out diseases which has less than 5 drugs
drug_shortlist = Filter(function(x) dim(x)[1] >= 1, drug_shortlist)
drug_shortlist_df = do.call(rbind, drug_shortlist)
drug_shortlist_df$Drug = tolower(drug_shortlist_df$Drug)
drug_shortlist_df$Drug = capitalize(drug_shortlist_df$Drug)
fwrite(drug_shortlist_df, file = file.path(dataFolder,"drug_shortlist.csv"))
save(drug_path,dis_path,drug_dis_path,drug_correlation,drug_shortlist,file=file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))
# load(file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))

drug_correlation_df = do.call(rbind, drug_correlation)
x = drug_correlation_df[drug_correlation_df$Correlation.Score >= 0.5 & drug_correlation_df$affectedPathway >= 80]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Scatter Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugCor = do.call(rbind, drug_correlation)
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
for (i in seq_along(drug_correlation)) {
    p[[i]] = ggplot(drug_correlation[[i]], aes(sample = Correlation.Score, col = Disease)) +
        stat_qq() +
        stat_qq_line() +
        labs(title = drug_correlation[[i]]$Disease) +
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

for (i in seq_along(drug_correlation)) {
    for (j in seq_along(drug_correlation[[i]])) {
        drug_correlation[[i]]$Dissimilarity.Score = NULL
        drug_correlation[[i]]$DrugPathway = NULL
        drug_correlation[[i]]$DiseasePathway = NULL
        drug_correlation[[i]]$Disease = NULL
        drug_correlation[[i]]$affectedPathway = NULL
    }
    names(drug_correlation[[i]])[[2]] = names(drug_correlation)[[i]]
}

# following drug(s) creates anomaly in the graph due to smaller density distribution in DrugGWAS data
# drug_correlation$`chronic obstructive pulmonary disease` = NULL
# 
# # following drug(s) creates anomaly in the graph due to smaller density distribution in DrugGWAS data
# # drug_correlation$`HIV infection` = NULL
# # drug_correlation$`celiac disease` = NULL
# # drug_correlation$`mucocutaneous lymph node syndrome` = NULL
# drug_correlation$`central nervous system cancer` = NULL
# drug_correlation$`Pick disease` = NULL
# drug_correlation$`Parkinson's disease` = NULL
# drug_correlation$`mucocutaneous lymph node syndrome` = NULL
# drug_correlation$`non-small cell lung carcinoma` = NULL
# drug_correlation$`celiac disease` = NULL
density.score = lapply(drug_correlation, melt)
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
# drug_correlation = do.call(rbind,drug_correlation) # list to df
drug_correlation = lapply(drug_correlation, function(x) na.omit(x))
drug_correlation = lapply(drug_correlation, data.frame) # need to convert back to data frame since following operations do not work on data table

# scaling with z score transformation
drug_cor_scaled = drug_correlation
for (i in seq_along(drug_correlation)) {
    drug_cor_scaled[[i]] = lapply(drug_correlation[[i]][3], function(x) scale(x, scale = TRUE))
    drug_cor_scaled[[i]]$durg = drug_correlation[[i]]$Drug
    drug_cor_scaled[[i]]$affectedPathway = drug_correlation[[i]]$affectedPathway
    #drug_cor_scaled[[i]]$disease = names(drug_cor_scaled)[[i]]
    drug_cor_scaled[[i]] = data.table(do.call(cbind, drug_cor_scaled[[i]]))
    drug_cor_scaled[[i]] = drug_cor_scaled[[i]][, c(2, 1, 3)]
    drug_cor_scaled[[i]]$V1 = as.numeric(drug_cor_scaled[[i]]$V1)
    names(drug_cor_scaled[[i]])[2] = names(drug_correlation[[i]])[2]
}

drug_cor_scaled$`celiac disease` = NULL #it gives spike on the graph due to less variabllity in the data
drug_cor_scaled$`non-small cell lung carcinoma` = NULL

drug_cor_scaled = lapply(drug_cor_scaled, melt)
for (i in seq_along(drug_cor_scaled)) {
    drug_cor_scaled[[i]]$variable = names(drug_cor_scaled)[[i]]
}
drug_cor_scaled = do.call(rbind, drug_cor_scaled)
drug_cor_scaled = data.table(drug_cor_scaled %>% drop_na())
drug_cor_scaled = drug_cor_scaled[, c(1, 3, 4, 2)]
names(drug_cor_scaled) = c("Drug", "Disease", "Correlation.Score", "affectedPathway")
drug_cor_scaled$Disease = as.character(drug_cor_scaled$Disease)
drug_cor_scaled$affectedPathway = as.numeric(drug_cor_scaled$affectedPathway)
drug_cor_scaled = drug_cor_scaled[order(drug_cor_scaled$Correlation.Score, - drug_cor_scaled$affectedPathway),]
drug_cor_scaled$Drug <- tolower(drug_cor_scaled$Drug)
drug_cor_scaled$Drug = capitalize(drug_cor_scaled$Drug)
#drug_cor_scaled$Disease = capitalize(drug_cor_scaled$Disease)


# jpeg(file=file.path(dataFolder,"densityPlots_DrugPdisease_CorrelationScore_scaled.jpeg"), width=2800, height=1980, res=200)
ggplot(drug_cor_scaled, aes(x = Correlation.Score, fill = Disease)) +
    geom_density(alpha = 0.25) +
    labs(title = "Distribution of Correlation Scores of all Drugs for each Diseases") +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 12)) +
    ylab("Density") +
    xlab("Correlation.Score")
dev.off()


# jpeg(file=file.path(dataFolder,"ScatterPlots_DrugPDisease_scaled_CorrelationScore.jpeg"), width=2800, height=1980, res=200)
ggplot(drug_cor_scaled, aes(x = affectedPathway, y = Correlation.Score, col = Disease)) +
    geom_point(size = 2, shape = 1) +
    labs(title = "Scatter Plots of Correlation Scores and affected pathways(%)") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Correlation Sore") +
    xlab("affected Pathway (%)")
dev.off()
