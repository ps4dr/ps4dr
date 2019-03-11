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

setwd("/home/memon/projects/msdrp/")

# load gwas.data data to get EFO IDs for diseases
load("./data/gwas.data.RData")
efo.ids = unique(gwas.data[, efo.id])

min.drugs.score = 0.2 #using 0.2 is a good trade-off to filter lower quality data points according to Open Targets
min.degs.score = 0.2 #using 0.2 is a good trade-off to filter lower quality data points according to Open Targets

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~Differentially Expressed Genes (DEGs)~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

DEGs = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
  tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter", 
            query = list(disease = efo.ids[i], datatype = "rna_expression", 
                         scorevalue_min = min.degs.score, size = 10000))
  tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
  if (length(tmp$data) > 0) {
    tmp = unique(as.data.table(flatten(tmp$data)))
    tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, ensembl.id = target.id, gene.symbol = target.gene_info.symbol, comparison = evidence.comparison_name, lfc = evidence.log2_fold_change.value, pval = evidence.resource_score.value)]
  }
}

DEGs = unique(DEGs)
length(unique(DEGs$efo.term))
length(unique(DEGs$gene.symbol))

save(DEGs, file="./data/DEGs.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~Drug-Disease relationships~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug2disease  = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter",
              query = list(disease = efo.ids[i], datatype = "known_drug", 
                           scorevalue_min = min.drugs.score, size = 10000))
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp = unique(as.data.table(flatten(tmp$data)))
        tmp[, chembl.id := sub("http://identifiers.org/chembl.compound/", "", drug.id)]
        tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, chembl.id, chembl.name = drug.molecule_name, chembl.type = drug.molecule_type, phase = drug.max_phase_for_all_diseases.numeric_index)]
    }
}

drug2disease = unique(drug2disease)
length(unique(drug2disease$efo.term))
length(unique(drug2disease$chembl.name))

save(drug2disease, file="./data/drug2disease.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~Therapeutic Area for each Diseases~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug2disease.therapeutic = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
  tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter", 
             query = list(disease = efo.ids[i], datatype = "rna_expression", 
                          scorevalue_min = min.degs.score, size = 10000))
  tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
  if (length(tmp$data) > 0) {
    tmp = unique(as.data.table(flatten(tmp$data)))
    tmp = tmp[, .(efo.id = disease.id,therapeutic.area = disease.efo_info.therapeutic_area.labels, efo.term = disease.efo_info.label)]
  }
}

save(drug2disease.therapeutic, file="./data/drug2disease.therapeutic.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


