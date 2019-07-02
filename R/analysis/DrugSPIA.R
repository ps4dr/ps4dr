#' 3rd script
#' summary:
#' This script calculates pathway enrichment analysis using Signaling Pathway Impact (SPIA) with:
#' 01: genes that are drug perturbed, differentially expressed and having SNPs from GWASs and/or
#' 02: genes that are drug perturbed and having SNPs from GWASs

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(tidyr)))

suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doParallel)))
suppressWarnings(suppressMessages(library(doSNOW)))
suppressWarnings(suppressMessages(library(doMC)))

no_cores <- detectCores() - 4

suppressWarnings(suppressMessages(library(SPIA)))

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####___Drug Data Preperation for SPIA___#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"geneID.RData"))
gene_id$entrezgene = gsub("^$", NA, gene_id$entrezgene)
gene_id = gene_id[which(! is.na(gene_id$entrezgene)),]
gene_id = data.table(gene_id)
gene_id = unique(gene_id[, c('entrezgene', 'ensembl_gene_id', 'hgnc_symbol')])
names(gene_id) = c("ENTREZ", "ensembl.id", "HGNC")
gene_id = gene_id[! duplicated(gene_id$ensembl.id),]
gene_id = gene_id[, c(1, 2)]

# get LINCS dataset
load(file.path(dataFolder,"L1000.RData"))
L1000 = L1000[, c(5, 1, 7)]
length(unique(L1000$chembl.id))
length(unique(L1000$ensembl.id))
L1000 = L1000[order(ensembl.id, decreasing = TRUE),]
L1000 = L1000[! duplicated(L1000[, c('ensembl.id', 'chembl.id')]),]
L1000 = merge(L1000, gene_id, by = "ensembl.id")

load(file.path(dataFolder,"drug2715details.RData"))
dmap = data.table(dmap[, c(1, 2, 3)])
dmap = dmap[phase == 4 | phase == 3 | phase == 2 | phase == 1] #673 Drugs
L1000 = merge(dmap[, c(1, 2)], L1000, by = "chembl.id")
# length(unique(L1000$chembl.id))

#' Users can select any of following two data sets 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~01: Drug Perturbed GWAS Genes~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load(file.path(dataFolder,"drugGWAS.genes50.padj1e-5.RData"))
# super_drugs = drugGWAS.genes[,c(2,3,4,9)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~02: Drug Perturbed Disease Genes (DEGs_GWAS)~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load(file.path(dataFolder,"drugPdisease_genes50.padj1e-5.RData"))
load(file.path(dataFolder,"drugPdisease_genes50.padj.RData"))
super_drugs = drugPdisease_genes[, c(2, 3, 4, 9)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' to select only those diseases from disease_genes50
load(file.path(dataFolder,"disease_genes50.RData"))
disease_genes = unique(disease_genes[, c(1, 2)])
super_drugs = merge(super_drugs, disease_genes, by.x = "efo.id", by.y = "efo.id.DEGs")
length(unique(super_drugs$efo.id))

# split multiple genes in same column to multiple rows
drug_genes = unique(super_drugs %>%
    mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>%
    unnest(ensembl.id))
#clean up special symbols
drug_genes$ensembl.id = gsub("\"", "", drug_genes$ensembl.id)
drug_genes$ensembl.id = gsub("c\\(", "", drug_genes$ensembl.id)
drug_genes$ensembl.id = gsub("\\)", "", drug_genes$ensembl.id)
drug_genes$ensembl.id = trimws(drug_genes$ensembl.id)

sprintf("Number of unique ChEMBL terms: %d", length(unique(drug_genes$chembl.name)))
sprintf("Number of unique EFO terms: %d", length(unique(drug_genes$efo.id)))
sprintf("Number of unique Genes: %d", length(unique(drug_genes$ensembl.id)))

#' Since we compare the Drug SPIA results with Disease SPIA
#' we will filter out all the entries from drug_genes which 
#' do not match the EFO.ID from Disease SPIA
#' 
load(file.path(dataFolder,"spia_output/spia_kegg_diseaseGenes50.RData"))
Disease47 = as.data.frame(names(spia_kegg))
names(Disease47) = "efo.term"
drug_genes = merge(Disease47, drug_genes, by.x = "efo.term", by.y = "efo.term.DEGs")

#' now merge with clinically trialled Drugs only from L1000 Data
drug_genes = merge(drug_genes, L1000, by = c('chembl.name', 'ensembl.id')) # merge with L1000
drug_genes = drug_genes[! duplicated(drug_genes[, c('chembl.id', 'ensembl.id', 'efo.term')]),]

#' count efo_chembl pair frequency, to remove rows frequencies less than 10
# tmp = count(drug_genes, c('efo.term','chembl.id'))
# tmp2 = subset(tmp, freq> 10)
# drug_genes = merge(tmp2,drug_genes, by="efo.term") # memory problem

topDiseases = split(drug_genes, drug_genes$efo.term) # split on diseases first
topDiseases_drug = lapply(topDiseases, function(x) split(x, x$chembl.name)) # split on drug inside each disease
topDiseases_drug = topDiseases_drug[lapply(topDiseases_drug, length) > 10] # remove diseases with very few drugs to test


### remove chemble.name, chemble.id and efo.term columns sicne they are not required any more
for (element in 1 : length(topDiseases_drug)) {
    for (subelem in 1 : length(topDiseases_drug[[element]])) {
        topDiseases_drug[[element]][[subelem]]$chembl.name = NULL
        topDiseases_drug[[element]][[subelem]]$chembl.id = NULL
        topDiseases_drug[[element]][[subelem]]$efo.id = NULL
        topDiseases_drug[[element]][[subelem]]$efo.term = NULL
        topDiseases_drug[[element]][[subelem]]$efo.term.y = NULL
    }
}

# save(topDiseases_drug,file = file.path(dataFolder,"drug_genes.list.drugGWAS.disease47.genes50.padj1e-10.RData"))
save(topDiseases_drug, file = file.path(dataFolder,"drug_genes.list.drugPdiseaseGenes50.padj.RData"))


####____Create named entity list for each drug_disease pairs___####

# load(file.path(dataFolder,"drug_genes.list.drugGWAS.disease47.genes50.padj1e-10.RData"))
load(file.path(dataFolder,"drug_genes.list.drugPdiseaseGenes50.padj.RData"))

lfc_drug_ensembl = topDiseases_drug
for (i in 1 : length(topDiseases_drug)) {
    for (j in 1 : length(topDiseases_drug[[i]])) {
        lfc_drug_ensembl[[i]][[j]] = setNames(as.numeric(topDiseases_drug[[i]][[j]][[2]]), as.character(topDiseases_drug[[i]][[j]][[1]]))
    }
}

lfc_drug_entrez = topDiseases_drug
for (i in 1 : length(topDiseases_drug)) {
    for (j in 1 : length(topDiseases_drug[[i]])) {
        lfc_drug_entrez[[i]][[j]] = setNames(as.numeric(topDiseases_drug[[i]][[j]][[2]]), as.character(topDiseases_drug[[i]][[j]][[3]]))
    }
}

lfc_drug_entrezID = topDiseases_drug
for (i in 1 : length(topDiseases_drug)) {
    for (j in 1 : length(topDiseases_drug[[i]])) {
        lfc_drug_entrezID[[i]][[j]] = setNames(as.numeric(topDiseases_drug[[i]][[j]][[2]]), as.character(gsub("^", "ENTREZID:", topDiseases_drug[[i]][[j]][[3]])))
    }
}

##_____Create a vector with all Gene universe to proxy Array Genes_________###
load(file.path(dataFolder,"geneID.RData"))
ensembl_all = unique(gene_id$ensembl.id)
entrez_all = unique(gene_id$ENTREZ)
entrezID_all = unique(gsub("^", "ENTREZID:", gene_id$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_drug_ensembl, ensembl_all, file = file.path(dataFolder,"lfc.drug.drugPdisease.diseaseGenes50.padj.ensembl.RData"))
save(lfc_drug_entrez, entrez_all, file = file.path(dataFolder,"lfc.drug.drugPdisease.diseaseGenes50.padj.entrez.RData"))
save(lfc_drug_entrezID, entrezID_all, file = file.path(dataFolder,"lfc.drug.drugPdisease.diseaseGenes50.padj.entrezID.RData"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~Drug SPIA~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~SPIA Function~~~~~~~~~~#
spia_fun = function(x){
    spia_drug_kegg = list()
    spia_drug_kegg = spia(de = x, all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}

#~~~~~~~~~Get Disease-Drugs list~~~~~~~~~#

# load(file.path(dataFolder,"lfc.drug.drugGWAS.diseaseGenes50.padj1e-5.entrez.RData"))
load(file.path(dataFolder,"lfc.drug.drugPdisease.diseaseGenes50.padj.entrez.RData"))
lfc_test = lfc_drug_entrez
#rm(lfc_drug_2)

##____Single Node______###

spia_kegg_drug = list()

for (diseases in 1 : length(lfc_test)) {
    spia_kegg_drug[[diseases]] = mclapply(lfc_test[[diseases]], spia_fun, mc.cores = no_cores)
    names(spia_kegg_drug)[[diseases]] = names(lfc_test)[[diseases]]
}

for (i in 1 : length(spia_kegg_drug)) {
    spia_kegg_drug[[i]] = Filter(function(x) ! is.null(x), spia_kegg_drug[[i]])
}

# save(spia_kegg_drug,file = file.path(dataFolder,"spia_output/spia_kegg_drug_Diseases_drugGWAS_nopar.RData"))
save(spia_kegg_drug, file = file.path(dataFolder,"spia_output/spia_kegg_drug_Diseases_drugPdisease_nopar.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
