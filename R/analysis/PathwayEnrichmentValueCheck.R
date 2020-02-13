#' additional script
#' summary:
#' This script calculates correation and dissimilarity scores for all drug-disease Gene pairs
#' This demonstrates how the results will differ if we use anti-correlation scores of Genes instead of disrupted pathways


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(RecordLinkage)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(Hmisc)))
suppressWarnings(suppressMessages(library(purrr)))
suppressWarnings(suppressMessages(library(tools)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(qqplotr)))

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~: load Disease Gene Expression Profiles :~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("/home/memon/projects/ps4dr/ps4dr/data/spia_input/lfc_disease_genes.RData") # without a filter for log-fold changes
load("/home/memon/projects/ps4dr/ps4dr/data/spia_input/lfc_2_disease_genes.RData") # log-fold changes >= |2|

# remove unnecessary columns
for (i in seq_along(lfc_efo)) {
  for (j in seq_along(lfc_efo[[i]])) {
    lfc_efo[[i]]$efo.id = NULL
    lfc_efo[[i]]$efo.term = NULL
    lfc_efo[[i]]$ENTREZ = NULL
    lfc_efo[[i]]$HGNC = NULL
  }
}

lfc_disease = lfc_efo

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#: load Drug perturbed Expression Profiles in Diseases :~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("/home/memon/projects/ps4dr/ps4dr/data/spia_input/lfc_drug_genes.RData")
lfc_drug = topDiseases_drug
rm(topDiseases_drug)


# remove unnecessary columns
for (element in 1 : length(lfc_drug)) {
  for (subelem in 1 : length(lfc_drug[[element]])) {
    lfc_drug[[element]][[subelem]]$ENTREZ = NULL
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#: remove diseases from lfc_rugs which are not in lfc_efo :#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Now, remove those diseases from lfc_efo which are not in drug_path
#' when lfc_drug smaller than lfc_disease

# lfc_disease = vector('list', length(lfc_drug)) # create list of lists
# names(lfc_disease) = names(lfc_drug)
# 
# for (i in 1 : length(lfc_efo)) {
#   if (names(lfc_efo)[[i]] %in% names(lfc_drug)) {
#     lfc_disease[[i]] = lfc_efo[[i]]
#     names(lfc_disease)[[i]] = names(lfc_efo)[[i]]
#   }
# }
# 
# lfc_disease = discard(lfc_disease, ~ all(is.na(.x)))
# rm(lfc_efo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#: remove diseases from lfc_efo which are not in lfc_drug :#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' when lfc_disease smaller than lfc_drug

lfc_drugs = vector('list', length(lfc_disease)) # create list of lists
names(lfc_drugs) = names(lfc_disease)

for (i in 1 : length(lfc_drug)) {
  if (names(lfc_drug)[[i]] %in% names(lfc_disease)) {
    lfc_drugs[[i]] = lfc_drug[[i]]
    names(lfc_drugs)[[i]] = names(lfc_drug)[[i]]
  }
}

lfc_drugs = discard(lfc_drugs, ~ all(is.na(.x)))
lfc_drug = lfc_drugs
lfc_disease$azoospermia = NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~: Calculate Correlation-Score :~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_dis_path = vector('list', length(lfc_drug)) # create list of lists
names(drug_dis_path) = names(lfc_drug)
drug_correlation = vector('list', length(lfc_drug)) # create list of lists
names(drug_correlation) = names(lfc_drug)

for (i in seq_along(lfc_drug)) {
  for (j in seq_along(lfc_drug[[i]])) {
    drug_dis_path[[i]][[j]] = merge(lfc_disease[[i]], lfc_drug[[i]][[j]], by = "ensembl.id")
    names(drug_dis_path[[i]])[[j]] = names(lfc_drug[[i]])[[j]]
    names(drug_dis_path[[i]][[j]]) = c("Genes", "Disease.Influence", "Drug.Influence")
    drug_dis_path[[i]][[j]]$Disease.Influence = ifelse(drug_dis_path[[i]][[j]]$Disease.Influence >= 0, 1, - 1)
    
    drug_correlation[[i]][[j]] = round(cor(drug_dis_path[[i]][[j]]$Drug.Influence, drug_dis_path[[i]][[j]]$Disease.Influence), 2)
    names(drug_correlation[[i]])[[j]] = names(drug_dis_path[[i]])[[j]]
    drug_correlation[[i]][[j]] = as.data.frame(drug_correlation[[i]][[j]])
    names(drug_correlation[[i]][[j]]) = "Correlation.Score"
    drug_correlation[[i]][[j]]$"Dissimilarity.Score" = round((sum(levenshteinDist(as.character(drug_dis_path[[i]][[j]]$Drug.Influence), as.character(drug_dis_path[[i]][[j]]$Disease.Influence))) * 100) / length(drug_dis_path[[i]][[j]]$Drug.Influence), 2)
    drug_correlation[[i]][[j]]$"DrugGenes" = length(drug_dis_path[[i]][[j]]$Drug.Influence)
    drug_correlation[[i]][[j]]$"DiseaseGenes" = length(lfc_disease[[i]]$ensembl.id)
    drug_correlation[[i]][[j]]$"affectedGenes" = round((drug_correlation[[i]][[j]]$"DrugGenes" / drug_correlation[[i]][[j]]$"DiseaseGenes") * 100, 2)
    drug_correlation[[i]][[j]]$"Disease" = names(lfc_disease)[[i]]
  }
}

#' create a single data.frame from all drugs in a disease

for (i in seq_along(drug_correlation)) {
  drug_correlation[[i]] = do.call(rbind, drug_correlation[[i]])
  drug_correlation[[i]]$Drug = rownames(drug_correlation[[i]])
  drug_correlation[[i]] = drug_correlation[[i]][, c(7, 6, 1, 2, 3, 4, 5)]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_correlation = do.call(rbind, drug_correlation) # list to df
drug_correlation = filter(drug_correlation, ! is.na(Correlation.Score)) # remove all rows with correlation score NA
drug_correlation = split(drug_correlation, drug_correlation$Disease) # df to list
drug_correlation = lapply(drug_correlation, data.table) # make all df to data table

#~~~~Remove any disease with only positive correlation scores for all drugs ~~~#
for (i in seq_along(drug_correlation)) {
  drug_correlation[[i]] = drug_correlation[[i]][! all(drug_correlation[[i]]$Correlation.Score >= 0)]
}

drug_correlation = Filter(function(x) dim(x)[1] >= 1, drug_correlation) # remove empty disease lists

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' filter out drugs from each disease with correlationscore greater than -0.1
drug_shortlist = drug_correlation
drug_shortlist = lapply(drug_shortlist, data.table)

#' remove drugs with less than the given thresholds
for (i in seq_along(drug_shortlist)) {
  drug_shortlist[[i]] = drug_shortlist[[i]][drug_shortlist[[i]]$Correlation.Score <= - 0.4 & drug_shortlist[[i]]$'affectedGenes' >= 50]
}

#' Order all lists in drug_shortlist
for (i in seq_along(drug_shortlist)) {
  drug_shortlist[[i]] = drug_shortlist[[i]][order(drug_shortlist[[i]][['Correlation.Score']],-drug_shortlist[[i]][['affectedGenes']]),]
}

#' filter out diseases which has less than 5 drugs
drug_shortlist = Filter(function(x) dim(x)[1] >= 1, drug_shortlist)
drug_shortlist_df = do.call(rbind, drug_shortlist)
drug_shortlist_df$Drug = tolower(drug_shortlist_df$Drug)
drug_shortlist_df$Drug = toTitleCase(drug_shortlist_df$Drug)
drug_shortlist_df$Disease = toTitleCase(drug_shortlist_df$Disease)
drug_shortlist_df = drug_shortlist_df[,c(2,1,3,4,5,6,7)]
names(drug_shortlist_df) = c("Disease","Drug","Correlation Score","Dissimilarity Score","Drug Influenced Genes","Disease Influenced Genes","Affected Genes(%)")

save(lfc_drug,lfc_disease,drug_dis_path,drug_correlation,drug_shortlist,file=file.path(dataFolder,"results/drugCorrelation_result_genes_all.RData"))
fwrite(drug_shortlist_df, file = file.path(dataFolder,"results/drug_shortlist_genes_all.csv"))

save(lfc_drug,lfc_disease,drug_dis_path,drug_correlation,drug_shortlist,file="/home/memon/projects/ps4dr/ps4dr/data/results/drugCorrelation_result_genes_lfc_2.RData")
fwrite(drug_shortlist_df, file = "/home/memon/projects/ps4dr/ps4dr/data/results/drug_shortlist_genes_lfc_2.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~: Scatter Plot :~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' remove few diseases which gives spike on the graph due to extremely less variabllity in the data
drug_correlation$`celiac disease` = NULL 
drug_correlation$`non-small cell lung carcinoma` = NULL
drug_correlation$`squamous cell carcinoma` = NULL

drugCor = do.call(rbind, drug_correlation)
drugCor = data.table(drugCor %>% drop_na())
drugCor = drugCor[order(drugCor$Correlation.Score, - drugCor$affectedGenes),]
drugCor$Drug <- tolower(drugCor$Drug)
drugCor$Drug = capitalize(drugCor$Drug)
drugCor$Disease = toTitleCase(drugCor$Disease)


jpeg(file = file.path(dataFolder, "results/figures/ScatterPlots_CorrelationScore_Genes.jpeg"), width = 3000, height = 1980, res = 200)
ggplot(drugCor, aes(x = affectedGenes, y = Correlation.Score, col = Disease)) +
  geom_point(size = 2, shape = 1) +
  labs(title = "Combined Scatter Plots of Drug's Correlation Scores and Affected Geness (%) in each Disease") +
  theme(legend.position = "bottom", legend.title = element_text(size = 10)) + theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Correlation Scores") +
  xlab("Affected Geness (%)")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~: Ferrero et al results :~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ferrero2 = read.csv("/home/memon/projects/ps4dr/ps4dr/data/13040_2018_171_MOESM2_ESM.csv")
ferrero3 = read.csv("/home/memon/projects/ps4dr/ps4dr/data/13040_2018_171_MOESM3_ESM.csv")

