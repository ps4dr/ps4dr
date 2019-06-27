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
gene.id$entrezgene = gsub("^$", NA, gene.id$entrezgene)
gene.id = gene.id[which(! is.na(gene.id$entrezgene)),]
gene.id = data.table(gene.id)
gene.id = unique(gene.id[, c('entrezgene', 'ensembl_gene_id', 'hgnc_symbol')])
names(gene.id) = c("ENTREZ", "ensembl.id", "HGNC")
gene.id = gene.id[! duplicated(gene.id$ensembl.id),]
gene.id = gene.id[, c(1, 2)]

# get LINCS dataset
load(file.path(dataFolder,"L1000.RData"))
L1000 = L1000[, c(5, 1, 7)]
length(unique(L1000$chembl.id))
length(unique(L1000$ensembl.id))
L1000 = L1000[order(ensembl.id, decreasing = TRUE),]
L1000 = L1000[! duplicated(L1000[, c('ensembl.id', 'chembl.id')]),]
L1000 = merge(L1000, gene.id, by = "ensembl.id")

load(file.path(dataFolder,"drug2715details.RData"))
dmap = data.table(dmap[, c(1, 2, 3)])
dmap = dmap[phase == 4 | phase == 3 | phase == 2 | phase == 1] #673 Drugs
L1000 = merge(dmap[, c(1, 2)], L1000, by = "chembl.id")
length(unique(L1000$chembl.id))

#' Users can select any of following two data sets 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~01: Drug Perturbed GWAS Genes~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load(file.path(dataFolder,"drugGWAS.genes50.padj1e-5.RData"))
# super.drugs = drugGWAS.genes[,c(2,3,4,9)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~02: Drug Perturbed Disease Genes (DEGs_GWAS)~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load(file.path(dataFolder,"drugPdisease.genes50.padj1e-5.RData"))
load(file.path(dataFolder,"drugPdisease.genes50.padj.RData"))
super.drugs = drugPdisease.genes[, c(2, 3, 4, 9)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' to select only those diseases from disease.genes50
load(file.path(dataFolder,"disease.genes50.RData"))
disease.genes = unique(disease.genes[, c(1, 2)])
super.drugs = merge(super.drugs, disease.genes, by.x = "efo.id", by.y = "efo.id.DEGs")
length(unique(super.drugs$efo.id))

# split multiple genes in same column to multiple rows
drug.genes = unique(super.drugs %>%
    mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>%
    unnest(ensembl.id))
#clean up special symbols
drug.genes$ensembl.id = gsub("\"", "", drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("c\\(", "", drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("\\)", "", drug.genes$ensembl.id)
drug.genes$ensembl.id = trimws(drug.genes$ensembl.id)

length(unique(drug.genes$chembl.name))
length(unique(drug.genes$efo.id))
length(unique(drug.genes$ensembl.id))

#' Since we compare the Drug SPIA results with Disease SPIA
#' we will filter out all the entries from drug.genes which 
#' do not match the EFO.ID from Disease SPIA
#' 
load(file.path(dataFolder,"spia_output/spia_kegg_disease47.genes50_results.RData"))
Disease47 = as.data.frame(names(spia_kegg))
names(Disease47) = "efo.term"
drug.genes = merge(Disease47, drug.genes, by.x = "efo.term", by.y = "efo.term.DEGs")

#' now merge with clinically trialled Drugs only from L1000 Data
drug.genes = merge(drug.genes, L1000, by = c('chembl.name', 'ensembl.id')) # merge with L1000
drug.genes = drug.genes[! duplicated(drug.genes[, c('chembl.id', 'ensembl.id', 'efo.term')]),]

#' count efo_chembl pair frequency, to remove rows frequencies less than 10
# tmp = count(drug.genes, c('efo.term','chembl.id'))
# tmp2 = subset(tmp, freq> 10)
# drug.genes = merge(tmp2,drug.genes, by="efo.term") # memory problem

topDiseases = split(drug.genes, drug.genes$efo.term) # split on diseases first
topDiseases.drug = lapply(topDiseases, function(x) split(x, x$chembl.name)) # split on drug inside each disease
topDiseases.drug = topDiseases.drug[lapply(topDiseases.drug, length) > 10] # remove diseases with very few drugs to test


### remove chemble.name, chemble.id and efo.term columns sicne they are not required any more
for (element in 1 : length(topDiseases.drug)) {
    for (subelem in 1 : length(topDiseases.drug[[element]])) {
        topDiseases.drug[[element]][[subelem]]$chembl.name = NULL
        topDiseases.drug[[element]][[subelem]]$chembl.id = NULL
        topDiseases.drug[[element]][[subelem]]$efo.id = NULL
        topDiseases.drug[[element]][[subelem]]$efo.term = NULL
        topDiseases.drug[[element]][[subelem]]$efo.term.y = NULL
    }
}

# save(topDiseases.drug,file = file.path(dataFolder,"drug.genes.list.drugGWAS.disease47.genes50.padj1e-10.RData"))
save(topDiseases.drug, file = file.path(dataFolder,"drug.genes.list.drugPdisease.genes50.padj.RData"))


####____Create named entity list for each drug_disease pairs___####

# load(file.path(dataFolder,"drug.genes.list.drugGWAS.disease47.genes50.padj1e-10.RData"))
load(file.path(dataFolder,"drug.genes.list.drugPdisease.genes50.padj.RData"))

lfc_drug_ensembl = topDiseases.drug
for (i in 1 : length(topDiseases.drug)) {
    for (j in 1 : length(topDiseases.drug[[i]])) {
        lfc_drug_ensembl[[i]][[j]] = setNames(as.numeric(topDiseases.drug[[i]][[j]][[2]]), as.character(topDiseases.drug[[i]][[j]][[1]]))
    }
}

lfc_drug_entrez = topDiseases.drug
for (i in 1 : length(topDiseases.drug)) {
    for (j in 1 : length(topDiseases.drug[[i]])) {
        lfc_drug_entrez[[i]][[j]] = setNames(as.numeric(topDiseases.drug[[i]][[j]][[2]]), as.character(topDiseases.drug[[i]][[j]][[3]]))
    }
}

lfc_drug_entrezID = topDiseases.drug
for (i in 1 : length(topDiseases.drug)) {
    for (j in 1 : length(topDiseases.drug[[i]])) {
        lfc_drug_entrezID[[i]][[j]] = setNames(as.numeric(topDiseases.drug[[i]][[j]][[2]]), as.character(gsub("^", "ENTREZID:", topDiseases.drug[[i]][[j]][[3]])))
    }
}

##_____Create a vector with all Gene universe to proxy Array Genes_________###
load(file.path(dataFolder,"geneID.RData"))
ensembl_all = unique(gene.id$ensembl.id)
entrez_all = unique(gene.id$ENTREZ)
entrezID_all = unique(gsub("^", "ENTREZID:", gene.id$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_drug_ensembl, ensembl_all, file = file.path(dataFolder,"lfc.drug.drugPdisease.disease47.genes50.padj.ensembl.RData"))
save(lfc_drug_entrez, entrez_all, file = file.path(dataFolder,"lfc.drug.drugPdisease.disease47.genes50.padj.entrez.RData"))
save(lfc_drug_entrezID, entrezID_all, file = file.path(dataFolder,"lfc.drug.drugPdisease.disease47.genes50.padj.entrezID.RData"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~Drug SPIA~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~SPIA Function~~~~~~~~~~#
spia_fun = function(x){
    spia_drug_kegg = list()
    spia_drug_kegg = spia(de = x, all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}

#~~~~~~~~~Get Disease-Drugs list~~~~~~~~~#

# load(file.path(dataFolder,"lfc.drug.drugGWAS.disease47.genes50.padj1e-5.entrez.RData"))
load(file.path(dataFolder,"lfc.drug.drugPdisease.disease47.genes50.padj.entrez.RData"))
lfc_test = lfc_drug_entrez
#rm(lfc_drug_2)

##____Single Node______###

spia_kegg_47D = list()

for (diseases in 1 : length(lfc_test)) {
    spia_kegg_47D[[diseases]] = mclapply(lfc_test[[diseases]], spia_fun, mc.cores = no_cores)
    names(spia_kegg_47D)[[diseases]] = names(lfc_test)[[diseases]]
}

for (i in 1 : length(spia_kegg_47D)) {
    spia_kegg_47D[[i]] = Filter(function(x) ! is.null(x), spia_kegg_47D[[i]])
}

# save(spia_kegg_47D,file = file.path(dataFolder,"spia_output/spia_kegg_47Diseases_drugGWAS_nopar.RData"))
save(spia_kegg_47D, file = file.path(dataFolder,"spia_output/spia_kegg_47Diseases_drugPdisease_nopar.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
