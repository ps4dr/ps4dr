#' additional script : PS4DR and Ferrero workflow comparison
#' summary:
#' This script calculates correation and dissimilarity scores for all drug-disease Gene pairs
#' This demonstrates how the results will differ if we use anti-correlation scores of Genes instead of disrupted pathways


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))

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
#~~~:  Comparison of PS4DR & Ferrero et al results  :~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
efo = fread("/home/memon/projects/ps4dr/ps4dr/data/EFO.csv")
efo = efo[, c(1, 2, 3)] # filter out unnecessary columns
efo$`Class ID` = gsub(".*efo/", "", efo$`Class ID`) # replace urls to get only EFO IDs
efo = efo[efo$`Class ID` %like% "EFO",] # filter out other disease IDs other than EFO
efo = unite(efo, Disease, c('Preferred Label', Synonyms), sep = "|", remove = FALSE)
efo = efo[, c(1, 2)]

# create new rows for each synonyms of EFO IDs 
efo = unique(efo %>%
               mutate(Disease = strsplit(as.character(Disease), "\\|")) %>%
               unnest(Disease))
efo$Disease = gsub("^$", NA, efo$Disease)
efo = efo[which(! is.na(efo$Disease)),]
efo$Disease = toTitleCase(efo$Disease)
names(efo) = c("efo.id","Disease")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load PS4DR results
load("/home/memon/projects/ps4dr/ps4dr/data/results/drugCorrelation_result.RData")

ps4dr_shortlist = do.call(rbind, drug_shortlist)
ps4dr_shortlist$Disease = toTitleCase(ps4dr_shortlist$Disease)

# ps4dr_allNegative = do.call(rbind, drug_correlation)
# ps4dr_allNegative = ps4dr_allNegative[Correlation.Score < 0 ]
# ps4dr_allNegative$Disease = toTitleCase(ps4dr_allNegative$Disease)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# map disease names to efo.id
tmp = unique(ps4dr_shortlist[,c(2)])
tmp = merge(tmp,efo,by="Disease")
tmp[2,2] = "EFO_0003869"
ps4dr_shortlist = merge(ps4dr_shortlist,tmp,by ="Disease")
ps4dr_shortlist = ps4dr_shortlist[,c(8,2,3,4,5,6,7)]

# tmp = unique(ps4dr_allNegative[,c(2)])
# tmp = merge(tmp,efo,by="Disease")
# tmp[6,2] = "EFO_0003869"
# tmp = unique(tmp)
# ps4dr_allNegative = merge(ps4dr_allNegative,tmp,by ="Disease")
# ps4dr_allNegative = ps4dr_allNegative[,c(8,2,3,4,5,6,7)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load Ferrero results
fer_true = fread("/home/memon/projects/ps4dr/ps4dr/data/13040_2018_171_MOESM2_ESM.csv")
fer_false = fread("/home/memon/projects/ps4dr/ps4dr/data/13040_2018_171_MOESM3_ESM.csv")
ferrero_all = rbind(fer_true,fer_false)
ferrero_all$efo.term = toTitleCase(ferrero_all$efo.term)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Comparison
comparison_shortlist = merge(ps4dr_shortlist,ferrero_all, by.x=c("efo.id","Drug"),by.y=c("efo.id","chembl.name"))
length(unique(comparison_shortlist$efo.id))
comparison_shortlist = split(comparison_shortlist, comparison_shortlist$efo.term)

# comparison_allNegative = merge(ps4dr_allNegative,ferrero_all,  by.x=c("efo.id","Drug"),by.y=c("efo.id","chembl.name"))
# length(unique(comparison_allNegative$efo.id))
# comparison_allNegative = split(comparison_allNegative, comparison_allNegative$efo.term)


save(comparison_shortlist,file = "/home/memon/projects/ps4dr/ps4dr/data/results/drug_comparison_ps4dr-ferrero.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Subset Melanoma Drug Prediction from Ferrero
fer_melanoma_true = fer_true[efo.term == "melanoma"]
fer_melanoma_false = fer_false[efo.term == "melanoma"]
fer_melanoma = rbind(fer_melanoma_true,fer_melanoma_false)

# Subset Melanoma Drug Prediction from PS4DR
ps4dr_melanoma_all = drug_correlation[["melanoma"]]
ps4dr_melanoma_all = ps4dr_melanoma_all[Correlation.Score < 0 ]
ps4dr_melanoma = drug_shortlist[["melanoma"]]

# Comparison Melanoma
melanoma_comparison_shortlist = merge(ps4dr_melanoma,fer_melanoma, by.x="Drug",by.y="chembl.name")
melanoma_comparison_allNegative = merge(ps4dr_melanoma_all,fer_melanoma, by.x="Drug",by.y="chembl.name")

save(melanoma_comparison_shortlist,melanoma_comparison_allNegative, file = "/home/memon/projects/ps4dr/ps4dr/data/results/_melanoma_drug_comparison_ps4dr-ferrero.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

