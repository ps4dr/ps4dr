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
suppressWarnings(suppressMessages(library(tools)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(qqplotr)))
suppressWarnings(suppressMessages(library(pROC)))

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
#~~~~~~~~~~~~~~~~: load Drug SPIA results :~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' load Kegg Drug SPIA file
load(file.path(dataFolder,"results/spia_output/spia_kegg_drug.RData"))
spia_drug = spia_kegg_drug
rm(spia_kegg_drug)

spia_drug = Filter(function(x) ! is.null(x), spia_drug) #delete empty df from list

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#

for (i in seq_along(spia_drug)) {
  spia_drug[[i]] = lapply(spia_drug[[i]], function(x) x[x$pNDE <= 0.05,])
  spia_drug[[i]] = spia_drug[[i]][lapply(spia_drug[[i]], length) > 1]
  spia_drug[[i]] = Filter(function(x) ! dim(x)[1] == 0, spia_drug[[i]])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_path_temp = vector('list', length(spia_drug)) # create list of lists
names(drug_path_temp) = names(spia_drug)
for (i in seq_along(spia_drug)) {
    for (j in seq_along(spia_drug[[i]])) {
        drug_path_temp[[i]][[j]] = spia_drug[[i]][[j]][, c(1, 11)] # use 2 for ID
        names(drug_path_temp[[i]])[[j]] = names(spia_drug[[i]])[[j]]
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~: load KEGG Disease SPIA results :~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' load Kegg Disease SPIA file
load(file.path(dataFolder,"results/spia_output/spia_kegg_disease.RData"))
spia_disease = spia_kegg
rm(spia_kegg)

#' Remove any disease pathway with p.value (pNDE) >= 0.05
spia_disease = lapply(spia_disease, function(x) x[x$pNDE <= 0.05,])
spia_disease = discard(spia_disease, ~ all(is.na(.x))) # remove empty diseases


#' Create disiase paths with only pathway name and their activity status
#' Remove any other diseases which are not in drug_path_temp
#' so, both drug_path_temp & dis_path are equivalent

dis_path = vector('list', length(drug_path_temp)) # create list of lists
names(dis_path) = names(drug_path_temp)

for (i in seq_along(drug_path_temp)) {
    for (j in seq_along(spia_disease)) {
        if (names(drug_path_temp)[[i]] == names(spia_disease)[j]) {
            dis_path[[i]] = spia_disease[[j]][, c(1, 11)]
        }
    }
}

#' remove empty disease frames
dis_path = discard(dis_path, ~ all(is.na(.x))) # remove empty diseases


#' Now, remove those diseases from drug_path_temp which are not in dis_path

drug_path = vector('list', length(dis_path)) # create list of lists
names(drug_path) = names(dis_path)

for (i in 1 : length(drug_path_temp)) {
  if (names(drug_path_temp)[[i]] %in% names(spia_disease)) {
    drug_path[[i]] = drug_path_temp[[i]]
    names(drug_path)[[i]] = names(drug_path_temp)[[i]]
  }
}

drug_path = discard(drug_path, ~ all(is.na(.x))) # remove empty diseases
rm(drug_path_temp,spia_disease,spia_drug)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~: Calculate Correlation-Score :~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
    drug_shortlist[[i]] = drug_shortlist[[i]][drug_shortlist[[i]]$Correlation.Score <= - 0.4 & drug_shortlist[[i]]$'affectedPathway' >= 50]
}

#' Order all lists in drug_shortlist
for (i in seq_along(drug_shortlist)) {
  drug_shortlist[[i]] = drug_shortlist[[i]][order(drug_shortlist[[i]][['Correlation.Score']],-drug_shortlist[[i]][['affectedPathway']]),]
}

#' filter out diseases which has less than 5 drugs
drug_shortlist = Filter(function(x) dim(x)[1] >= 1, drug_shortlist)
drug_shortlist_df = do.call(rbind, drug_shortlist)
drug_shortlist_df$Drug = tolower(drug_shortlist_df$Drug)
drug_shortlist_df$Drug = toTitleCase(drug_shortlist_df$Drug)
drug_shortlist_df$Disease = toTitleCase(drug_shortlist_df$Disease)
drug_shortlist_df = drug_shortlist_df[,c(2,1,3,4,5,6,7)]
names(drug_shortlist_df) = c("Disease","Drug","Correlation Score","Dissimilarity Score","Drug Influenced Pathway","Disease Influenced Pathway","Affected Pathway(%)")

save(drug_path,dis_path,drug_dis_path,drug_correlation,drug_shortlist,file=file.path(dataFolder,"results/drugCorrelation_result.RData"))
fwrite(drug_shortlist_df, file = file.path(dataFolder,"results/drug_shortlist.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~: Cancer Drugs Ratio :~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' load drug details (phase=I,II,III,IV) from ChEMBL (28/01/2020)
tmp = read.csv(file.path(dataFolder,"chembl_1234.csv", sep = "\t"))
tmp = tmp[,c(1,2,3)]

#' load all 673 use case drugs
load(file.path(dataFolder,"usecase_drugs.RData"))
dmap = dmap[,c(1,2)]

tmp2 = merge(dmap,tmp,by.x="chembl.id", by.y="ChEMBL.ID")
tmp2 = tmp2[,c(3,4)]
chembl673 = unique(tmp2 %>% # split multiple genes in same column to multiple rows
                     mutate(Synonyms = strsplit(as.character(Synonyms), "\\|")) %>%
                     unnest(Synonyms))
chembl673$Synonyms = trimws(chembl673$Synonyms)
chembl673$Synonyms = toupper(chembl673$Synonyms)
chembl673 = unique(chembl673)
chembl673X=chembl673[,c(2,1)]
names(chembl673X) = c("Name","Synonyms")
tmp = rbind(chembl673,chembl673X)
dmap673 = unique(as.data.table(tmp[,c(1)]))

#' load all Cancer Drugs from National Cancer Institute (28/01/2020)
cancerDrugs = read.csv(file.path(dataFolder, "cancerdrugs_NCI_28012020.txt"))
cancerDrugs = unique(cancerDrugs)

cancerDrugs = unique(cancerDrugs %>% # split multiple genes in same column to multiple rows
                       mutate(cancerdrug = strsplit(as.character(cancerdrug), "\\(")) %>%
                       unnest(cancerdrug))
cancerDrugs$cancerdrug = gsub("\\)", "", cancerDrugs$cancerdrug)
cancerDrugs$cancerdrug = trimws(cancerDrugs$cancerdrug)
cancerDrugs$cancerdrug = toupper(cancerDrugs$cancerdrug)
cancerDrugs = unique(cancerDrugs)


#' merge cancer drugs with all 673 use case drugs (and their synonyms)
alldrugs_diseaseMap = merge(dmap673,cancerDrugs,by.x = "V1",by.y="cancerdrug")
sprintf("Cancer Drug Ratio in the all Use case Drugs: %.2f", length(unique(alldrugs_diseaseMap$V1))/673*100)

#' merge canser drugs with all 115 shortlisted drugs
topdrugs = unique(drug_shortlist_df[,c(1)])
topdrugs_diseaseMap = merge(topdrugs,alldrugs_diseaseMap, by.x="Drug",by.y = "V1")
sprintf("Cancer Drug Ratio in the Top Drug list: %.2f", length(unique(topdrugs_diseaseMap$Drug))/115*100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~: Scatter Plot :~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' remove few diseases which gives spike on the graph due to extremely less variabllity in the data
drug_correlation$`celiac disease` = NULL
drug_correlation$`non-small cell lung carcinoma` = NULL
drug_correlation$`squamous cell carcinoma` = NULL

drugCor = do.call(rbind, drug_correlation)
drugCor = data.table(drugCor %>% drop_na())
drugCor = drugCor[order(drugCor$Correlation.Score, - drugCor$affectedPathway),]
drugCor$Drug <- tolower(drugCor$Drug)
drugCor$Drug = capitalize(drugCor$Drug)
drugCor$Disease = toTitleCase(drugCor$Disease)


jpeg(file = file.path(dataFolder, "results/figures/ScatterPlots_CorrelationScore.jpeg"), width = 3000, height = 1980, res = 200)
ggplot(drugCor, aes(x = affectedPathway, y = Correlation.Score, col = Disease)) +
    geom_point(size = 2, shape = 1) +
    labs(title = "Combined Scatter Plots of Drug's Correlation Scores and Affected Pathways (%) in each Disease") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) + theme(legend.title=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Correlation Scores") +
    xlab("Affected Pathways (%)")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~: interactive scatterplot :~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

splot <- ggplot(drugCor, aes(x = affectedPathway, y = Correlation.Score, col = Disease,label=Drug)) +
  geom_point(size = 2, shape = 1) +
  labs(title = "Combined Scatter Plots of Drug's Correlation Scores and Affected Pathways (%) in each Disease") +
  theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +theme(legend.title=element_blank()) +
  ylab("Correlation Scores") +
  xlab("Affected Pathways (%)")
splot_int = ggplotly(splot,tooltip =c("Disease","Drug","Correlation.Score"))
htmlwidgets::saveWidget(as_widget(splot_int), file.path(dataFolder, "results/figures/ScatterPlots_CorrelationScore.html"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~: Density Plot :~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load(file.path(dataFolder,"drugCorrelation_result.RData"))
#
# for (i in seq_along(drug_correlation)) {
#     for (j in seq_along(drug_correlation[[i]])) {
#         drug_correlation[[i]]$Dissimilarity.Score = NULL
#         drug_correlation[[i]]$DrugPathway = NULL
#         drug_correlation[[i]]$DiseasePathway = NULL
#         drug_correlation[[i]]$Disease = NULL
#         drug_correlation[[i]]$affectedPathway = NULL
#     }
#     names(drug_correlation[[i]])[[2]] = names(drug_correlation)[[i]]
# }
#
# density.score = lapply(drug_correlation, melt)
# density.score = do.call(rbind, density.score)
# names(density.score) = c("Drug", "Diseases", "Correlation_Score")
#
# # jpeg(file=file.path(dataFolder, "results/figures/densityPlots_allDiseases.jpeg"), width=2800, height=1980, res=200)
# ggplot(density.score, aes(x = Correlation_Score, fill = Diseases)) +
#     geom_density(alpha = 0.25) +
#     labs(title = "Distribution of Correlation Coefficient Scores of all Drugs for each Diseases") +
#     theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
#     theme(plot.title = element_text(hjust = 0.5)) +theme(legend.title=element_blank()) +
#     ylab("Density") +
#     xlab("Correlation Coefficient Scores")
# dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~: Density Plot With Standization :~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
load(file.path(dataFolder,"results/drugCorrelation_result.RData"))

#' Standization with scale function (Z-score normalization)
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

drug_cor_scaled$`celiac disease` = NULL #gives spike on the graph due to less variabllity in the data
drug_cor_scaled$`non-small cell lung carcinoma` = NULL
#squamous cell carcinoma get excluded itself since it has only 1 value for all drugs

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
drug_cor_scaled$Disease = tools::toTitleCase(drug_cor_scaled$Disease)

jpeg(file=file.path(dataFolder, "results/figures/densityPlots_CorrelationScore_scaled.jpeg"), width=3000, height=1980, res=190)
ggplot(drug_cor_scaled, aes(x = Correlation.Score, fill = Disease)) +
    geom_density(alpha = 0.25) +
    labs(title = "Distribution of Z-score Normalized Correlation Scores of all Drugs for each Diseases") +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
    theme(legend.title=element_blank()) +
    ylab("Density") +
    xlab("Normalized Correlation Scores")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~: Q-Q Plot With Standization :~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

jpeg(file=file.path(dataFolder, "results/figures/Q-Q_Plots_CorrelationScore_scaled.jpeg"), width=3000, height=1980, res=190)
ggplot(data = drug_cor_scaled, mapping = aes(sample = Correlation.Score, color = Disease, fill = Disease)) +
  stat_qq_band(alpha=0.5) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~ Disease) + facet_grid(rows = 5) +
  theme(legend.position = "none", legend.title = element_text(size = 10)) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~ Density Plots for use case scenario Diseases ~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

DisUseCase = c("Alzheimer's Disease","Breast Carcinoma","Melanoma","Pancreatic Carcinoma")
case = as.data.frame(DisUseCase)
case = merge(drug_cor_scaled,case, by.x="Disease", by.y="DisUseCase")

# jpeg(file=file.path(dataFolder, "results/figures/densityPlots_CorrelationScore_scaled_useCase.jpeg"), width=3000, height=1980, res=200)
p1 <- ggplot(case, aes(x = Correlation.Score, fill = Disease)) +
  labs(title = "Distribution of Z-score normalized correlation scores of drugs\nfor 4 use case scenario diseases") +
  geom_density(alpha = 0.25) +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
  theme(legend.title=element_blank()) + ylab("Density") +
  xlab("Normalized Correlation Scores")
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~ Density Plots for 4 Diseases have no drugs in shortlist ~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# DisnoShortlst = c("Amyotrophic Lateral Sclerosis","Asthma","Atopic Eczema","Autism","Chronic Lymphocytic Leukemia","Hypertriglyceridemia","Mucocutaneous Lymph Node Syndrome","Psoriasis","Ulcerative Colitis")
DisnoShortlst = c("Asthma","Mucocutaneous Lymph Node Syndrome","Psoriasis","Ulcerative Colitis")
xcase = as.data.frame(DisnoShortlst)
xcase = merge(drug_cor_scaled,xcase, by.x="Disease", by.y="DisnoShortlst")

# jpeg(file=file.path(dataFolder, "results/figures/densityPlots_CorrelationScore_scaled_noshortlist.jpeg"), width=3000, height=1980, res=200)
p2 <- ggplot(xcase, aes(x = Correlation.Score, fill = Disease)) +
  labs(title = "Distribution of Z-score normalized correlation scores of drugs\nfor 4 diseases which have no shortlisted drugs") +
  geom_density(alpha = 0.25) +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
  theme(legend.title=element_blank()) + ylab("Density") +
  xlab("Normalized Correlation Scores")
# dev.off()

jpeg(file=file.path(dataFolder, "results/figures/densityPlots_CorrelationScore_scaled_usecaseVsnoshortlist.jpeg"), width=3200, height=1280, res=180)
plot_grid(p1,p2)
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# twoDiseases = c("Alzheimer's Disease","Psoriasis")
twoDiseases = c("Melanoma","Psoriasis")
casexcase = as.data.frame(twoDiseases)
casexcase = merge(drug_cor_scaled,casexcase, by.x="Disease", by.y="twoDiseases")

jpeg(file=file.path(dataFolder, "results/figures/densityPlots_CorrelationScore_scaled_Melanoma_Psoriasis.jpeg"), width=3000, height=1980, res=200)
ggplot(casexcase, aes(x = Correlation.Score, fill = Disease)) +
  geom_density(alpha = 0.25) +
  labs(title = "Distribution of Z-score Normalized Correlation Scores of all Drugs for Melanoma and Psoriasis") +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10)) +
  theme(legend.title=element_blank()) + ylab("Density") +
  xlab("Normalized Correlation Scores")
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~ ROC Curve for predicted drugs ~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

predictedDrugs = fread(file.path(dataFolder, "results/drug_shortlist.csv"))
predictedDrugs = predictedDrugs[,c(2,1,3)]
names(predictedDrugs) = c("chembl.name","efo.term", "correlationScore")
predictedDrugs$chembl.name = toupper(predictedDrugs$chembl.name)
predictedDrugs$efo.term = tolower(predictedDrugs$efo.term)
predictedDrugs$efo.term = capitalize(predictedDrugs$efo.term)

# load drugPurturbed genesets list to map predictedDrugs with exisiting indications
load(file.path(dataFolder,"results/drugPdisease_genes.RData"))

drugPdisease_genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
drugPdisease_genes$efo.term = tolower(drugPdisease_genes$efo.term)
drugPdisease_genes$efo.term = capitalize(drugPdisease_genes$efo.term)

predictedDrugs = merge(predictedDrugs,drugPdisease_genes, by=c("chembl.name","efo.term"))

# Even though we have used OpenTarget for mapping with clinical trials info, we did manual curation for all
# our predicted drug-disease pairs by checking them in the  clinicaltirals.gov
# git link for manually curated clinical trial info: https://github.com/ps4dr/ps4dr/blob/master/data/Predicted-Drugs-with-clinical-Trials-Info.csv
# Below, we change 6 additional drug-disease association with the updated information from the above table

predictedDrugs[103,13] = TRUE; predictedDrugs[5,13] = TRUE; predictedDrugs[3,13] = TRUE; predictedDrugs[76,13] = TRUE; predictedDrugs[29,13] = TRUE; predictedDrugs[102,13] = TRUE

predictedDrugs = predictedDrugs[ order(correlationScore),]
predictedDrugs$correlationScore = predictedDrugs$correlationScore*-1
scale01 = function(x){(x-min(x)) / (max(x)-min(x))}

predictedDrugs$correlationScoreSclaed = scale01(predictedDrugs$correlationScore)
predictedDrugs = predictedDrugs[,c(1,2,3,15,13)]


preds <- predictedDrugs[, -log10(correlationScoreSclaed)]
labls <- as.numeric(predictedDrugs[, existing.indication])

drugROC <- roc(response = labls, predictor = preds, algorithm = 2, ci = TRUE, ci.method = "bootstrap", boot.n = 1000, parallel = TRUE, progress = "none")
print(drugROC)
save(ROCresult, file=file.path(dataFolder,"results/ROC_result.RDATA"))

ci.specificity <- ci.sp(drugROC, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
ci.sensivity <- ci.se(drugROC, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")

jpeg(file=file.path(dataFolder, "results/figures/predictedDrugs_ROC.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
par(pty = "s")
plot(smooth(drugROC), main = paste("AUC =", round(drugROC$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
plot(ci.sensivity, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(ci.specificity, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(smooth(drugROC), add = TRUE, col = "#ff6600", lwd = 3)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

