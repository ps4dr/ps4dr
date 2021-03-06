#' additional script for CREEDS data cleaning
#' summary:
#' This script clean the CREEDS data fop PS4DR pipeline

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

load(file.path(dataFolder,"geneID_v97.RData"))
gene_id$ENTREZ = gsub("^$", NA, gene_id$ENTREZ)
gene_id = gene_id[which(! is.na(gene_id$ENTREZ)),]

#~~~~~~~~~~~: CREEDS Disease data: Melanoma :~~~~~~~~~~~~~#

creeds_melanoma = fread(file.path(dataFolder, "creeds_input/melanoma_disease_data.tsv",sep = "\t"))
creeds_melanoma = creeds_melanoma[,c(2,4,5)]
names(creeds_melanoma) = c("efo.id","HGNC","lfc")
creeds_melanoma$efo.id = "EFO_0000756"
lfc_comb = creeds_melanoma
lfc_comb = merge(lfc_comb, gene_id, by = "HGNC")
lfc_comb = lfc_comb[! duplicated(lfc_comb[, c('efo.id', 'ENTREZ')]),]
lfc_comb = lfc_comb[,c(5,2,3,4,1)]
names(lfc_comb)
save(lfc_comb, file = file.path(dataFolder,"results/CREEDs/lfc_disease_creeds.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~: CREEDS Drug data: Bortezomib :~~~~~~~~~~~~~#

creeds_bortezomib = fread(file.path(dataFolder, "creeds_input/bortezomib_drug_data.tsv",sep = "\t"))
creeds_bortezomib = creeds_bortezomib[,c(2,4,5)]
names(creeds_bortezomib) = c("chembl.id","HGNC","lfc")
creeds_bortezomib$chembl.id = "CHEMBL325041"
creeds_bortezomib = merge(creeds_bortezomib,gene_id,by="HGNC")
creeds_bortezomib = creeds_bortezomib[,c(2,5,3)]
L1000 = creeds_bortezomib
save(L1000, file = file.path(dataFolder,"results/CREEDs/lfc_drug_creeds.RData"))

save(spia_kegg_creeds_melanoma, spia_kegg_drug_creeds, file = file.path(dataFolder, "results/CREEDs/spia_kegg_creeds.RData"))
save(drug_path,dis_path,drug_dis_path,drug_correlation, file = file.path(dataFolder, "results/CREEDs/drugCorrelation_result_creeds.RData"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
