#' 5th script
#' summary:
#' create pathway activities with the combination of two drugs
#' then use 2drugs combination pathways to calculate correlation score 

library(data.table)
library(tidyr)
library(RecordLinkage)
library(ggplot2)
library(purrr)

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

#' load results from pervious script (4-CheckDistribution.R)
# load(file.path(dataFolder,"drugCorraltion.drugGWAS.RData"))
load(file.path(dataFolder,"drugCorraltion.drugPdisease.RData"))

#' filter out drugs from each disease with correlationscore greater than 0
drug.shortlist = drug.Correlation
drug.shortlist = lapply(drug.shortlist, data.table)

for (i in seq_along(drug.shortlist)) {
    drug.shortlist[[i]] = drug.shortlist[[i]][drug.shortlist[[i]]$Correlation.Score < 0]
}

#' filter out diseases which has less than 1 drugs
drug.shortlist = Filter(function(x) dim(x)[1] >= 1, drug.shortlist)

#' create separate data frame for each drug
for (i in seq_along(drug.shortlist)) {
    drug.shortlist[[i]] = split(drug.shortlist[[i]], drug.shortlist[[i]]$Drug)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~load KEGG Drug SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load(file.path(dataFolder,"spia_output/spia_kegg_47Diseases_drugGWAS_nopar.RData"))
load(file.path(dataFolder,"spia_output/spia_kegg_47Diseases_drugPdisease_nopar.RData"))

spia_drug_kegg = spia_kegg_47D
rm(spia_kegg_47D)
spia_drug_kegg = Filter(function(x) ! is.null(x), spia_drug_kegg) #delete empty df from list

# remove empty drug frames
for (i in 1 : length(spia_drug_kegg)) {
    spia_drug_kegg[[i]] = spia_drug_kegg[[i]][lapply(spia_drug_kegg[[i]], length) > 1]
}

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#
for (i in seq_along(spia_drug_kegg)) {
    spia_drug_kegg[[i]] = lapply(spia_drug_kegg[[i]], function(x) x[x$pNDE <= 0.05,])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' since combination of drug pair is a computationally expensive
#' we will test it for the drugs that appeared in our drug.shortlist
#' So, fisrt we filter diseases from "spia_drug_kegg" which are not in drug.shortlist
#' Then, we filter the drugs in those idseases that are not in drug.shortlist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 1. Filter Diseases from "spia_drug_kegg" which are not in drug.shortlist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

spia_drug_kegg_short = vector('list', length(spia_drug_kegg))
names(spia_drug_kegg_short) = names(spia_drug_kegg)

for (i in 1 : length(spia_drug_kegg)) {
  if (names(spia_drug_kegg)[[i]] %in% names(drug.shortlist)) {
    spia_drug_kegg_short[[i]] = spia_drug_kegg[[i]]
    names(spia_drug_kegg_short)[[i]] = names(spia_drug_kegg)[[i]]
  }
}

#' remove empty disease frames
spia_drug_kegg_short = discard(spia_drug_kegg_short, ~ all(is.na(.x)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 2. Filter Drugs from Diseases in "spia_drug_kegg" which are not in drug.shortlist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugNegCor = vector('list', length(drug.shortlist)) # create list of lists
names(drugNegCor) = names(drug.shortlist)

for (i in 1 : length(spia_drug_kegg_short)) {
    for (j in 1 : length(spia_drug_kegg_short[[i]])) {
        if (names(spia_drug_kegg_short[[i]])[[j]] %in% names(drug.shortlist[[i]])) {
            drugNegCor[[i]][[j]] = spia_drug_kegg_short[[i]][[j]]
        }
    }
    drugNegCor[[i]] = Filter(function(x) ! is.null(x), drugNegCor[[i]]) # filter out empty drug table
    names(drugNegCor[[i]]) = names(drug.shortlist[[i]]) # name drugs
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Create data frame with afftected pathways and affect direction by the drugs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug.path = vector('list', length(drugNegCor)) # create list of lists
names(drug.path) = names(drugNegCor)

for (i in seq_along(drugNegCor)) {
    for (j in seq_along(drugNegCor[[i]])) {
        drug.path[[i]][[j]] = drugNegCor[[i]][[j]][, c(1, 11)] # use 2 for ID, 1 for Name
        names(drug.path[[i]])[[j]] = names(drugNegCor[[i]])[[j]]
    }
}

#' filter drug path with only 1 disease, since they can't be used for combination
drug.path = Filter(function(x) length(x) > 1, drug.path)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Create drug combinations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

comb1 = list() # combination of data
comb2 = list() # combination of drug names
set.seed(123)
for (i in seq_along(drug.path)) {
    comb1[[i]] = combn(drug.path[[i]], 2)
    comb2[[i]] = combn(names(drug.path[[i]]), 2)
}

rm(drug.dis.path, drug.shortlist, drug.Correlation,spia_drug_kegg,spia_drug_kegg_short)

drug.comb.path = vector('list', length(drug.path)) # create list of lists
names(drug.comb.path) = names(drug.path)

rm(drug.path)

for (i in 1 : length(comb1)) {
    for (j in 1 : (length(comb1[[i]]) / 2)) {

        drug.comb.path[[i]][[j]] = data.table(merge(comb1[[i]][, j][1], comb1[[i]][, j][2], by = "Name", all = TRUE)) #Name/ID
        names(drug.comb.path[[i]])[[j]] = paste(comb2[[i]][, j][1], comb2[[i]][, j][2], sep = "_")
        # create a dummy column by copmaring disease and drug affects columns
        drug.comb.path[[i]][[j]][, Status.update :=
        ifelse(is.na(drug.comb.path[[i]][[j]]$Status.x), drug.comb.path[[i]][[j]]$Status.y,
        ifelse(is.na(drug.comb.path[[i]][[j]]$Status.y), drug.comb.path[[i]][[j]]$Status.x,
        ifelse(drug.comb.path[[i]][[j]]$Status.x != drug.comb.path[[i]][[j]]$Status.y, "Neutral", drug.comb.path[[i]][[j]]$Status.y)))]
        drug.comb.path[[i]][[j]]$Status.x = NULL
        drug.comb.path[[i]][[j]]$Status.y = NULL
    }
}

# save(drug.comb.path,file=file.path(dataFolder,"drug.comb.path.drugGWAS.RData"))
save(drug.comb.path, file = file.path(dataFolder,"drug.comb.path.drugPdisease.RData"))
rm(comb1, comb2)
load(file.path(dataFolder,"drug.comb.path.drugPdisease.RData"))
# load(file.path(dataFolder,"drug.comb.path.drugGWAS.RData"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' discard diseases from dis.path which are not in drug.comb.path
a = setdiff(names(dis.path), names(drug.comb.path))
for (i in seq_along(dis.path)) {
    if (names(dis.path)[i] %in% a) {
        dis.path[[i]] = NULL
    }
}

drug.dis.path = vector('list', length(drug.comb.path)) # create list of lists
names(drug.dis.path) = names(drug.comb.path)
drug.Correlation = vector('list', length(drug.comb.path)) # create list of lists
names(drug.Correlation) = names(drug.comb.path)

for (i in seq_along(drug.comb.path)) {
    for (j in seq_along(drug.comb.path[[i]])) {
        drug.dis.path[[i]][[j]] = merge(dis.path[[i]], drug.comb.path[[i]][[j]], by = "Name")
        names(drug.dis.path[[i]])[[j]] = names(drug.comb.path[[i]])[[j]]
        names(drug.dis.path[[i]][[j]]) = c("Pathways", "Disease.Influence", "Drug.Influence")
        drug.dis.path[[i]][[j]]$Disease.Influence = ifelse(drug.dis.path[[i]][[j]]$Disease.Influence == "Activated", 1, - 1)
        drug.dis.path[[i]][[j]]$Drug.Influence = ifelse(drug.dis.path[[i]][[j]]$Drug.Influence == "Activated", 1, - 1)


        drug.Correlation[[i]][[j]] = cor(drug.dis.path[[i]][[j]]$Drug.Influence, drug.dis.path[[i]][[j]]$Disease.Influence)
        names(drug.Correlation[[i]])[[j]] = names(drug.dis.path[[i]])[[j]]
        drug.Correlation[[i]][[j]] = as.data.frame(drug.Correlation[[i]][[j]])
        names(drug.Correlation[[i]][[j]]) = "Correlation.Score"
        drug.Correlation[[i]][[j]]$"Dissimilarity.Score" = (sum(levenshteinDist(as.character(drug.dis.path[[i]][[j]]$Drug.Influence), as.character(drug.dis.path[[i]][[j]]$Disease.Influence))) * 100) / length(drug.dis.path[[i]][[j]]$Drug.Influence)
        drug.Correlation[[i]][[j]]$"DrugPathway" = length(drug.dis.path[[i]][[j]]$Drug.Influence)
        drug.Correlation[[i]][[j]]$"DiseasePathway" = length(dis.path[[i]]$Name) #change for ID/Name
        drug.Correlation[[i]][[j]]$"affectedPathway" = round((drug.Correlation[[i]][[j]]$"DrugPathway" / drug.Correlation[[i]][[j]]$"DiseasePathway") * 100, 2)
        drug.Correlation[[i]][[j]]$"Disease" = names(dis.path)[[i]]
    }
}

##~~~~~create a single data.frame from all drugs in a disease~~~##

#x=do.call(rbind,drug.Correlation[["Alzheimer's disease"]])
for (i in seq_along(drug.Correlation)) {
    drug.Correlation[[i]] = do.call(rbind, drug.Correlation[[i]])
    drug.Correlation[[i]]$Drug = rownames(drug.Correlation[[i]])
    drug.Correlation[[i]] = drug.Correlation[[i]][, c(7, 6, 1, 2, 3, 4, 5)]
}

save(drug.dis.path, drug.Correlation, file = file.path(dataFolder,"drugCombination.correlationScore.drugPdisease.RData"))

load(file.path(dataFolder,"drugCombination.correlationScore.drugPdisease.RData"))
x = drug.Correlation[["lung carcinoma"]][order(drug.Correlation[["lung carcinoma"]]$Correlation.Score, - drug.Correlation[["lung carcinoma"]]$affectedPathway),]
x = drug.Correlation[[5]][order(- drug.Correlation[[5]]$affectedPathway, drug.Correlation[[5]]$Correlation.Score),]
x = data.table(x)
x$Correlation.Score = round(x$Correlation.Score, 2)
x = x[x$Correlation.Score <= - 0.5 & x$affectedPathway >= 80]
x = x[order(- x$affectedPathway, x$Correlation.Score),]
x = x[, c(1, 3, 7)]
fwrite(x, file = file.path(dataFolder,"drugComb.breastCancer.csv"))
x$Drug = tolower(x$Drug)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Scatter Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugCor = do.call(rbind, drug.Correlation)
drugCor$Drug = as.factor(drugCor$Drug)
drugCor = data.table(drugCor %>% drop_na())
#' remove few extreme correlation scores
# drugCor = drugCor[Correlation.Score > -1 & Correlation.Score < 1]

# jpeg(file=file.path(dataFolder,"ScatterPlots_combination_Diseases_DrugGWAS.jpeg", width=2800, height=1980, res=200))
jpeg(file = file.path(dataFolder,"ScatterPlots_combination_Diseases_DrugPdisease.jpeg"), width = 2800, height = 1980, res = 200)
ggplot(drugCor, aes(x = affectedPathway, y = Correlation.Score, col = Disease)) +
    geom_point(size = 2, shape = 1) +
    labs(title = "Scatter Plots of Correlation Scores and affected pathways(%)") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Correlation Sore") +
    xlab("affected Pathway (%)")
dev.off()
# save(drugCor, file = file.path(dataFolder,"drugCombcorScoreDrugGWAS.RData"))
length(unique(drugCor$Disease))
save(drugCor, file = file.path(dataFolder,"drugCombcorScoreDrugPdisease.RData"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
