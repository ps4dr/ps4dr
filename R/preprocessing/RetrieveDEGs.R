#' 4th script
#' summary:
#' 01: Retrieve Gene expression data from the Open Targets by using their API service
#'     for every diseases with the EFO IDs from gwas.data data
#' 02: Retrieve Drug-Disease relationships and clinical trail phase inofrmations
#' 03: Retrieve Therapeutic area mapping for all diseases as well
#' last accessed on 11 March, 2019

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(httr)))
suppressWarnings(suppressMessages(library(jsonlite)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doParallel)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))

registerDoParallel(parallel::detectCores() - 1)

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

# load gwas.data data to get EFO IDs for diseases
load(file.path(dataFolder,"GWASs.RData"))
efo.ids = unique(GWASs[, efo.id])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~Retrieve Differentially Expressed Genes (DEGs)~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

DEGs = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter",
    query = list(disease = efo.ids[i], datatype = "rna_expression",
    scorevalue_min = - Inf, size = 10000))
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp = unique(as.data.table(flatten(tmp$data)))
        tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, ensembl.id = target.id, gene.symbol = target.gene_info.symbol, comparison = evidence.comparison_name, lfc = evidence.log2_fold_change.value, pval = evidence.resource_score.value)]
    }
    try(print(paste("log of", i, "=", log(i))))
}


DEGs = unique(DEGs)
sprintf("Number of unique EFO terms: %d", length(unique(DEGs$efo.term)))
sprintf("Number of unique Gene Symbol: %d", length(unique(DEGs$gene.symbol)))

save(DEGs, file = file.path(dataFolder, "DEGs.RData"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~Retrieve Drug-Disease relationships~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug2disease = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter",
    query = list(disease = efo.ids[i], datatype = "known_drug",
    scorevalue_min = - Inf, size = 10000))
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp = unique(as.data.table(flatten(tmp$data)))
        tmp[, chembl.id := sub("http://identifiers.org/chembl.compound/", "", drug.id)]
        tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, chembl.id, chembl.name = drug.molecule_name, chembl.type = drug.molecule_type, phase = drug.max_phase_for_all_diseases.numeric_index)]
    }
}

drug2disease = unique(drug2disease)
sprintf("Number of unique EFO terms: %d", length(unique(drug2disease$efo.term)))
sprintf("Number of unique ChEMBL Name: %d", length(unique(drug2disease$chembl.name)))

save(drug2disease, file = file.path(dataFolder, "drug2disease.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


