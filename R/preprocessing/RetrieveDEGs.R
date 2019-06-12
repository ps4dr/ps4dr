#' 4th script
#' summary:
#' 01: Retrieve Gene expression data from the Open Targets by using their API service
#'     for every diseases with the EFO IDs from gwas.data data
#' 02: Retrieve Drug-Disease relationships and clinical trail phase inofrmations
#' 03: Retrieve Therapeutic area mapping for all diseases as well
#' last accessed on 11 March, 2019

library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

#####################################################################
#TODO: Change to the directory where you cloned this repository
setwd("/home/memon/projects/msdrp/")
#####################################################################

# load gwas.data data to get EFO IDs for diseases
load("./data/gwas.data.RData")
efo.ids = unique(gwas.data[, efo.id])

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
}

DEGs = unique(DEGs)
length(unique(DEGs$efo.term))
length(unique(DEGs$gene.symbol))

save(DEGs, file = "./data/DEGs.RData")


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
        tmp[, chembl.id :  = sub("http://identifiers.org/chembl.compound/", "", drug.id)]
        tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, chembl.id, chembl.name = drug.molecule_name, chembl.type = drug.molecule_type, phase = drug.max_phase_for_all_diseases.numeric_index)]
    }
}

drug2disease = unique(drug2disease)
length(unique(drug2disease$efo.term))
length(unique(drug2disease$chembl.name))

save(drug2disease, file = "./data/drug2disease.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~Retrieve Therapeutic Area for each Diseases~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug2disease.therapeutic = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter",
    query = list(disease = efo.ids[i], datatype = "rna_expression",
    scorevalue_min = - Inf, size = 10000))
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp = unique(as.data.table(flatten(tmp$data)))
        tmp = tmp[, .(efo.id = disease.id, therapeutic.area = disease.efo_info.therapeutic_area.labels, efo.term = disease.efo_info.label)]
    }
}

drug2disease.therapeutic = unique(as.data.frame(drug2disease.therapeutic))

# create new rows for each synonyms of EFO IDs 
drug2disease.therapeutic = unique(drug2disease.therapeutic %>%
    mutate(therapeutic.area = strsplit(as.character(therapeutic.area), ",")) %>%
    unnest(therapeutic.area))
# cleaning up symbols
drug2disease.therapeutic$therapeutic.area = gsub("^c\\(", "", drug2disease.therapeutic$therapeutic.area)
drug2disease.therapeutic$therapeutic.area = gsub("\\)", "", drug2disease.therapeutic$therapeutic.area)
drug2disease.therapeutic$therapeutic.area = gsub("\"", "", drug2disease.therapeutic$therapeutic.area)

length(unique(drug2disease.therapeutic$efo.term))
length(unique(drug2disease.therapeutic$therapeutic.area))

save(drug2disease.therapeutic, file = "./data/drug2disease.therapeutic.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


