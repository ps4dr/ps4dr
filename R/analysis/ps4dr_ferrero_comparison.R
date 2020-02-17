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
#ps4dr_shortlist = ps4dr_shortlist[,c(8,2,3,4,5,6,7)]

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
comparison_shortlist_list = split(comparison_shortlist, comparison_shortlist$efo.id)

# comparison_allNegative = merge(ps4dr_allNegative,ferrero_all,  by.x=c("efo.id","Drug"),by.y=c("efo.id","chembl.name"))
# length(unique(comparison_allNegative$efo.id))
# comparison_allNegative = split(comparison_allNegative, comparison_allNegative$efo.term)

ps4dr_shortlist = split(ps4dr_shortlist,ps4dr_shortlist$Disease)

comparison_table <- data.frame(matrix(ncol = 4, nrow = 12))
names(comparison_table) = c("Disease","# PS4DR_Drug","# match_with_Ferrero","PS4DR_Ferrero_Match(%)")

for (i in 1 : length(comparison_shortlist_list)) {
  for (j in 1 : length(ps4dr_shortlist)) {  
    #if (names(comparison_shortlist_list)[[i]] == names(ps4dr_shortlist)[[j]]) {
    if (names(comparison_shortlist_list)[[i]] == ps4dr_shortlist[[j]][["efo.id"]][[1]]) {
      comparison_table[[1]][[i]] = names(ps4dr_shortlist)[[j]]
      comparison_table[[2]][[i]] = nrow(ps4dr_shortlist[[j]])
      comparison_table[[3]][[i]] = nrow(comparison_shortlist_list[[i]])
      comparison_table[[4]][[i]] = round(nrow(comparison_shortlist_list[[i]])/nrow(ps4dr_shortlist[[j]])*100,2)
    }
  }
}

save(comparison_shortlist,file = "/home/memon/projects/ps4dr/ps4dr/data/results/drug_comparison_ps4dr-ferrero.RData")
write.csv(comparison_table,file = "/home/memon/projects/ps4dr/ps4dr/data/results/drug_comparison_table_ps4dr-ferrero.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

