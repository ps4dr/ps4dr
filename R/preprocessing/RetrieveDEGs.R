#' 2nd script
#' summary:
#' 01: Retrieve Gene expression data from the Open Targets by using their API service
#'     for every diseases with the EFO IDs from gwas.data data
#' 02: Retrieve Drug-Disease relationships and clinical trail phase inofrmations
#' last accessed on 11 March, 2019

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(httr)))
suppressWarnings(suppressMessages(library(jsonlite)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))

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
efo_ids = unique(GWASs[, efo.id])
rm(GWASs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~Retrieve Differentially Expressed Genes (DEGs)~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!file.exists(file.path(dataFolder, "DEGs.RData"))){
  DEGs = data.table()
  tryCatch(for(i in 1:length(efo_ids)){
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter", 
              query = list(disease = efo_ids[i], datatype = "rna_expression"  ,scorevalue_min = -Inf, size = 10000))
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
      cat(sprintf("retrieving differential gene expression data for: %s, index #%d \n", efo_ids[i], i))
      tmp2 = data.table(efo.id = tmp$data$disease$id, efo.term = tmp$data$disease$efo_info$label, ensembl.id = tmp$data$target$id,
                        gene.symbol = tmp$data$target$gene_info$symbol, comparison = tmp$data$evidence$comparison_name, 
                        lfc = tmp$data$evidence$log2_fold_change$value, pval = tmp$data$evidence$resource_score$value)
      tmp2 = unique(tmp2)
      DEGs = rbind(DEGs,tmp2)
    }
  },
  error=function(e) 1)
  rm(tmp,tmp2)
  DEGs = unique(DEGs)
  cat(sprintf("Number of unique EFO terms: %d\n", length(unique(DEGs$efo.term))))
  cat(sprintf("Number of unique Gene Symbol: %d\n", length(unique(DEGs$gene.symbol))))
  
  save(DEGs, file = file.path(dataFolder, "DEGs.RData"))
} else {cat(sprintf("~~ DEGs file already exists, not downloading again. ~~\n"))}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~Retrieve Drug-Disease relationships~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!file.exists(file.path(dataFolder, "drug2disease.RData"))){
  drug2disease = data.table()
  tryCatch(for(i in 1:length(efo_ids)){
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter",
              query = list(disease = efo_ids[i], datatype = "known_drug",
                           scorevalue_min = - Inf, size = 10000))
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
      cat(sprintf("retrieving drug information for : %s, index #%d \n", efo_ids[i], i ))
      tmp2 = data.table(efo.id = tmp$data$disease$id, efo.term = tmp$data$disease$efo_info$label, chembl.id = tmp$data$drug$id,
                        chembl.name = tmp$data$drug$molecule_name, chembl.type = tmp$data$drug$molecule_type, 
                        phase = tmp$data$drug$max_phase_for_all_diseases$numeric_index)
      tmp2[, chembl.id := sub("http://identifiers.org/chembl.compound/", "", chembl.id)]
      tmp2 = unique(tmp2)
      drug2disease = rbind(drug2disease,tmp2)
    }
  },
  error=function(e) 1)
  rm(tmp,tmp2)
  drug2disease = unique(drug2disease)
  cat(sprintf("Number of unique EFO terms: %d\n", length(unique(drug2disease$efo.term))))
  cat(sprintf("Number of unique ChEMBL Name: %d\n", length(unique(drug2disease$chembl.name))))
  
  save(drug2disease, file = file.path(dataFolder, "drug2disease.RData"))
} else {cat(sprintf("~~ drug2disease file already exists, not downloading again. ~~\n"))}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


