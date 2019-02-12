#' 4th script
#' summary:
#' Retrieve the Open Targets Gene expression data for every diseases by using the EFO IDs from STOPGAP data
#' Retrieve Therapeutic area mapping for all diseases as well

library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

# options
# permissive
min.drugs.score = -Inf
min.degs.score = -Inf
## strict
#min.drugs.score = 0.5
#min.degs.score = 0.2

setwd("/home/memon/projects/msdrp/")

# read STOPGAP data to get EFO IDs for diseases

stopgap = fread("./data/stopgap.tsv")
#stopgap2 = fread("./data/stopgap_snp.tsv")

# diseases to process
efo.ids = unique(stopgap[, efo.id])


## drugs
# loop thorugh diseases
opentargets.drugs  = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {

    # query API
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "known_drug", scorevalue_min = min.drugs.score, size = 10000))

    # extract and check data
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp = unique(as.data.table(flatten(tmp$data)))
        tmp[, chembl.id := sub("http://identifiers.org/chembl.compound/", "", drug.id)]
        tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, chembl.id, chembl.name = drug.molecule_name, chembl.type = drug.molecule_type, phase = drug.max_phase_for_all_diseases.numeric_index)]
    }
}

# export
fwrite(opentargets.drugs, "./data/opentargets.drugs.tsv", sep = "\t")


##### Differentially Expressed Genes (DEGs) ######

# loop thorugh diseases
opentargets.degs = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
    
    # query API
    tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter", query = list(disease = efo.ids[i], datatype = "rna_expression", scorevalue_min = min.degs.score, size = 10000))

    # extract and check data
    tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
    if (length(tmp$data) > 0) {
        tmp = unique(as.data.table(flatten(tmp$data)))
        tmp = tmp[, .(efo.id = disease.id, efo.term = disease.efo_info.label, ensembl.id = target.id, gene.symbol = target.gene_info.symbol, comparison = evidence.comparison_name, lfc = evidence.log2_fold_change.value, pval = evidence.resource_score.value)]
    }
}

# export
fwrite(opentargets.degs, "./data/opentargets.degs.tsv", sep = "\t")

length(unique(opentargets.degs$efo.term))
length(unique(opentargets.degs$gene.symbol))

length(unique(opentargets.drugs$efo.term))
length(unique(opentargets.drugs$chembl.name))


##### Therapeutic Area for each Diseases ######
# loop thorugh diseases
opentargets.degs.therapeutic = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
  tmp = GET("https://api.opentargets.io/v3/platform/public/evidence/filter", 
             query = list(disease = efo.ids[i], datatype = "rna_expression", 
                          scorevalue_min = min.degs.score, size = 10000))
  
  tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
  if (length(tmp$data) > 0) {
    tmp = unique(as.data.table(flatten(tmp$data)))
    tmp = tmp[, .(efo.id = disease.id,therapeutic.area = disease.efo_info.therapeutic_area.labels, efo.term = disease.efo_info.label)]
  }
}

fwrite(opentargets.degs.therapeutic, "./data/opentargets.degs.therapeutic.tsv", sep = "\t")



harmonizome.ctd = foreach(i = seq(efo.ids), .combine = rbind) %dopar% {
  
  # query API
  tmp = GET("Harmonizome/api/1.0/gene", query = list(disease = efo.ids[i], datatype = "gene"))
  
  # extract and check data
  tmp = fromJSON(content(tmp, type = "text", encoding = "UTF-8"))
  if (length(tmp$data) > 0) {
    tmp = unique(as.data.table(flatten(tmp$data)))
  }
}



