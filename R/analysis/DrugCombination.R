#' 5th script
#' summary:
#' create pathway activities with the combination of two drugs
#' then use 2drugs combination pathways to calculate correlation score 

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(RecordLinkage)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(purrr)))
suppressWarnings(suppressMessages(library(stringr)))

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
load(file.path(dataFolder,"results/drugCorrelation_result.RData"))

#' filter out drugs from each disease with correlationscore greater than 0
drug_shortlist = drug_correlation
drug_shortlist = lapply(drug_shortlist, data.table)

for (i in seq_along(drug_shortlist)) {
    drug_shortlist[[i]] = drug_shortlist[[i]][drug_shortlist[[i]]$Correlation.Score < 0]
}

#' filter out diseases which has less than 1 drugs
drug_shortlist = Filter(function(x) dim(x)[1] >= 1, drug_shortlist)

#' create separate data frame for each drug
for (i in seq_along(drug_shortlist)) {
    drug_shortlist[[i]] = split(drug_shortlist[[i]], drug_shortlist[[i]]$Drug)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~load KEGG Drug SPIA results~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"results/spia_output/spia_kegg_drug.RData"))

spia_kegg_drug = Filter(function(x) ! is.null(x), spia_kegg_drug) #delete empty df from list

# remove empty drug frames
for (i in 1 : length(spia_kegg_drug)) {
    spia_kegg_drug[[i]] = spia_kegg_drug[[i]][lapply(spia_kegg_drug[[i]], length) > 1]
}

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#
for (i in seq_along(spia_kegg_drug)) {
    spia_kegg_drug[[i]] = lapply(spia_kegg_drug[[i]], function(x) x[x$pNDE <= 0.05,])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' since combination of drug pair is a computationally expensive
#' we will test it for the drugs that appeared in our drug_shortlist
#' So, fisrt we filter diseases from "spia_kegg_drug" which are not in drug_shortlist
#' Then, we filter the drugs in those idseases that are not in drug_shortlist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 1. Filter Diseases from "spia_kegg_drug" which are not in drug_shortlist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

spia_kegg_drug_short = vector('list', length(spia_kegg_drug))
names(spia_kegg_drug_short) = names(spia_kegg_drug)

for (i in 1 : length(spia_kegg_drug)) {
  if (names(spia_kegg_drug)[[i]] %in% names(drug_shortlist)) {
    spia_kegg_drug_short[[i]] = spia_kegg_drug[[i]]
    names(spia_kegg_drug_short)[[i]] = names(spia_kegg_drug)[[i]]
  }
}

#' remove empty disease frames
spia_kegg_drug_short = discard(spia_kegg_drug_short, ~ all(is.na(.x)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' 2. Filter Drugs from Diseases in "spia_kegg_drug" which are not in drug_shortlist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugNegCor = vector('list', length(drug_shortlist)) # create list of lists
names(drugNegCor) = names(drug_shortlist)

for (i in 1 : length(spia_kegg_drug_short)) {
    for (j in 1 : length(spia_kegg_drug_short[[i]])) {
        if (names(spia_kegg_drug_short[[i]])[[j]] %in% names(drug_shortlist[[i]])) {
            drugNegCor[[i]][[j]] = spia_kegg_drug_short[[i]][[j]]
        }
    }
    drugNegCor[[i]] = Filter(function(x) ! is.null(x), drugNegCor[[i]]) # filter out empty drug table
    names(drugNegCor[[i]]) = names(drug_shortlist[[i]]) # name drugs
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Create data frame with afftected pathways and affect direction by the drugs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_path = vector('list', length(drugNegCor)) # create list of lists
names(drug_path) = names(drugNegCor)

for (i in seq_along(drugNegCor)) {
    for (j in seq_along(drugNegCor[[i]])) {
        drug_path[[i]][[j]] = drugNegCor[[i]][[j]][, c(1, 11)] # use 2 for ID, 1 for Name
        names(drug_path[[i]])[[j]] = names(drugNegCor[[i]])[[j]]
    }
}

#' filter drug path with only 1 disease, since they can't be used for combination
drug_path = Filter(function(x) length(x) > 1, drug_path)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Create drug combinations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

comb1 = list() # combination of data
comb2 = list() # combination of drug names
set.seed(123)
for (i in seq_along(drug_path)) {
    comb1[[i]] = combn(drug_path[[i]], 2)
    comb2[[i]] = combn(names(drug_path[[i]]), 2)
}

rm(drug_dis_path, drug_shortlist, drug_correlation,spia_kegg_drug,spia_kegg_drug_short)

drug_comb_path = vector('list', length(drug_path)) # create list of lists
names(drug_comb_path) = names(drug_path)

rm(drug_path)

for (i in 1 : length(comb1)) {
    for (j in 1 : (length(comb1[[i]]) / 2)) {

        drug_comb_path[[i]][[j]] = data.table(merge(comb1[[i]][, j][1], comb1[[i]][, j][2], by = "Name", all = TRUE)) #Name/ID
        names(drug_comb_path[[i]])[[j]] = paste(comb2[[i]][, j][1], comb2[[i]][, j][2], sep = "_")
        # create a dummy column by copmaring disease and drug affects columns
        drug_comb_path[[i]][[j]][, Status.update :=
        ifelse(is.na(drug_comb_path[[i]][[j]]$Status.x), drug_comb_path[[i]][[j]]$Status.y,
        ifelse(is.na(drug_comb_path[[i]][[j]]$Status.y), drug_comb_path[[i]][[j]]$Status.x,
        ifelse(drug_comb_path[[i]][[j]]$Status.x != drug_comb_path[[i]][[j]]$Status.y, "Neutral", drug_comb_path[[i]][[j]]$Status.y)))]
        drug_comb_path[[i]][[j]]$Status.x = NULL
        drug_comb_path[[i]][[j]]$Status.y = NULL
    }
}

# save(drug_comb_path,file=file.path(dataFolder,"drug_comb_path.drugGWAS.RData"))
save(drug_comb_path, file = file.path(dataFolder,"results/drug_comb_path_drugPdisease.RData"))
rm(comb1, comb2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' discard diseases from dis_path which are not in drug_comb_path

load(file.path(dataFolder,"results/drug_comb_path_drugPdisease.RData"))

a = setdiff(names(dis_path), names(drug_comb_path))
for (i in seq_along(dis_path)) {
    if (names(dis_path)[i] %in% a) {
        dis_path[[i]] = NULL
    }
}

drug_dis_path = vector('list', length(drug_comb_path)) # create list of lists
names(drug_dis_path) = names(drug_comb_path)
drug_correlation = vector('list', length(drug_comb_path)) # create list of lists
names(drug_correlation) = names(drug_comb_path)

for (i in seq_along(drug_comb_path)) {
    for (j in seq_along(drug_comb_path[[i]])) {
        drug_dis_path[[i]][[j]] = merge(dis_path[[i]], drug_comb_path[[i]][[j]], by = "Name")
        names(drug_dis_path[[i]])[[j]] = names(drug_comb_path[[i]])[[j]]
        names(drug_dis_path[[i]][[j]]) = c("Pathways", "Disease.Influence", "Drug.Influence")
        drug_dis_path[[i]][[j]]$Disease.Influence = ifelse(drug_dis_path[[i]][[j]]$Disease.Influence == "Activated", 1, - 1)
        drug_dis_path[[i]][[j]]$Drug.Influence = ifelse(drug_dis_path[[i]][[j]]$Drug.Influence == "Activated", 1, - 1)


        drug_correlation[[i]][[j]] = cor(drug_dis_path[[i]][[j]]$Drug.Influence, drug_dis_path[[i]][[j]]$Disease.Influence)
        names(drug_correlation[[i]])[[j]] = names(drug_dis_path[[i]])[[j]]
        drug_correlation[[i]][[j]] = as.data.frame(drug_correlation[[i]][[j]])
        names(drug_correlation[[i]][[j]]) = "Correlation.Score"
        drug_correlation[[i]][[j]]$"Dissimilarity.Score" = (sum(levenshteinDist(as.character(drug_dis_path[[i]][[j]]$Drug.Influence), as.character(drug_dis_path[[i]][[j]]$Disease.Influence))) * 100) / length(drug_dis_path[[i]][[j]]$Drug.Influence)
        drug_correlation[[i]][[j]]$"DrugPathway" = length(drug_dis_path[[i]][[j]]$Drug.Influence)
        drug_correlation[[i]][[j]]$"DiseasePathway" = length(dis_path[[i]]$Name) #change for ID/Name
        drug_correlation[[i]][[j]]$"affectedPathway" = round((drug_correlation[[i]][[j]]$"DrugPathway" / drug_correlation[[i]][[j]]$"DiseasePathway") * 100, 2)
        drug_correlation[[i]][[j]]$"Disease" = names(dis_path)[[i]]
    }
}

##~~~~~create a single data.frame from all drugs in a disease~~~##

#x=do.call(rbind,drug_correlation[["Alzheimer's disease"]])
for (i in seq_along(drug_correlation)) {
    drug_correlation[[i]] = do.call(rbind, drug_correlation[[i]])
    drug_correlation[[i]]$Drug = rownames(drug_correlation[[i]])
    drug_correlation[[i]] = drug_correlation[[i]][, c(7, 6, 1, 2, 3, 4, 5)]
}

save(drug_dis_path, drug_correlation, file = file.path(dataFolder,"results/drugCombination_result_drugPdisease.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##~~~~~~~~~~~~~~~~~~~~~Scatter Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"results/drugCombination_result_drugPdisease.RData"))

drugCor = do.call(rbind, drug_correlation)
drugCor$Drug = as.factor(drugCor$Drug)
drugCor = data.table(drugCor %>% drop_na())
#' remove few extreme correlation scores
# drugCor = drugCor[Correlation.Score > -1 & Correlation.Score < 1]

jpeg(file = file.path(dataFolder, "results/figures/ScatterPlots_combination_Diseases.jpeg"), width = 2800, height = 1980, res = 200)
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
save(drugCor, file = file.path(dataFolder,"results/drugCombination_CorrelationScore.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drugCor_Shortlist = drugCor[Correlation.Score <= -.5 & affectedPathway >= 80 ]
drugCor_breastCancer_1 = drugCor[Disease == "breast carcinoma" & Correlation.Score <= -.4 & affectedPathway >= 50 ] # 464 Drug Combinations
drugCor_breastCancer_2 = drugCor[Disease == "breast carcinoma" & Correlation.Score <= -.5 & affectedPathway >= 80 ] # 26 Drug Combinations
drug_comb_bc= cbind(str_split_fixed(drugCor_breastCancer_2$Drug, "_", 2),drugCor_breastCancer_2)
drug_comb_bc = drug_comb_bc[,c(1,2,5,6,7,8,9)]
names(drug_comb_bc) = c("Drug_1","Drug_2","Correlation Score","Dissimilarity Score","Drug Pathway","Disease Pathway","Affected Pathway(%)")
drug_comb_bc$Drug_1 = tolower(drug_comb_bc$Drug_1)
drug_comb_bc$Drug_1 = toTitleCase(drug_comb_bc$Drug_1)
drug_comb_bc$Drug_2 = tolower(drug_comb_bc$Drug_2)
drug_comb_bc$Drug_2 = toTitleCase(drug_comb_bc$Drug_2)


fwrite(drug_comb_bc, file = file.path(dataFolder,"results/drugCombination_shortlist_breastCancer.csv"))

# fwrite(drug_comb_bc, file = "/home/memon/projects/ps4dr/ps4dr/data/results/drugCombination_shortlist_breastCancer.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
