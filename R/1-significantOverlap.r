#' 1st script
#' summary:
#' This script calculates significant overlap between/among different data sets using Fisher's exact tests.
#' Three different overlaps are calculated in this script (Figure 1. "Gene set intersection")
#' 01: GWAS & DEGs data :> disease.genes,
#' 02: Drug Perturbed Genes & DEGs & GWAS data :> drugPdisease.genes,
#' 03: Drug Perturbed Genes & GWAS data :> drugGWAS.genes.

library(EnsDb.Hsapiens.v86)

library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


#####################################################################
#TODO: Change to the directory where you cloned this repository
setwd("/home/memon/projects/msdrp/")
#####################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~Enrichment Calculation Function for Two Gene Sets using Fisher's Exact Test~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SignificantOverlap <- function(geneSet1, geneSet2, universe) {
  geneSet1 = unique(geneSet1)
  geneSet2 = unique(geneSet2)
  x = list(intersect(geneSet1,geneSet2))
  a = sum(geneSet1 %in% geneSet2)
  b = length(geneSet1) - a
  c = length(geneSet2) - a
  d = length(universe) - a - b - c
  fisher = fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  data.table(DEGs = length(geneSet1), GWASs = length(geneSet2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~Datasets~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("./data/GWASs.RData")
tmp = dplyr::count(GWASs, efo.id)
tmp2 = subset(tmp, n > 50)
GWASs = merge(GWASs,tmp2,by="efo.id")

load("./data/DEGs.RData")
tmp = dplyr::count(DEGs, efo.id)
tmp2 = subset(tmp, n > 50)
DEGs = merge(DEGs,tmp2,by="efo.id")

#################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____DEGs to GWAS Genes >>> DISEASE Genes_____####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.

# create gene universes for Fisher's test
ensembl.ids = unique(keys(EnsDb.Hsapiens.v86))

# split data table by disease
GWASs.list = split(GWASs, GWASs$efo.id)
DEGs.list = split(DEGs, DEGs$efo.id)

## Signifcant overlap calculation
disease.genes <- foreach (i = seq(DEGs.list), .combine = rbind, .errorhandling = "remove") %do% {
  efo.id.DEGs = names(DEGs.list)[i]
  foreach (j = seq(GWASs.list), .combine = rbind, .errorhandling = "remove") %dopar% {
    efo.id.GWASs = names(GWASs.list)[j]
    tmp = SignificantOverlap(DEGs.list[[efo.id.DEGs]]$ensembl.id, GWASs.list[[efo.id.GWASs]]$ensembl.id, ensembl.ids)
    tmp = cbind(efo.id.DEGs, efo.id.GWASs, tmp)
    tmp = merge(tmp, unique(GWASs[, .(efo.id, efo.term)]), by.x = "efo.id.GWASs", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
    tmp = merge(tmp, unique(DEGs[, .(efo.id, efo.term)]), by.x = "efo.id.DEGs", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".GWASs", ".DEGs"))
    setcolorder(tmp, c("efo.id.DEGs", "efo.term.DEGs", "efo.id.GWASs", "efo.term.GWASs", "DEGs", "GWASs", "overlap", "universe","commonGenes", "odds.ratio", "p.value"))
  }
}

disease.genes = unique(disease.genes)
# correct p-values
disease.genes[p.value == 0, p.value := 3e-324]
disease.genes = disease.genes[order(p.value), ]
disease.genes = disease.genes[overlap > 0]
disease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
disease.genes = disease.genes[p.adjusted <= 0.05]
# add dummy variable
disease.genes[, same.disease := ifelse(efo.id.DEGs == efo.id.GWASs, TRUE, FALSE)]
save(disease.genes,file="./data/disease.genes50.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##___________Drugs to DISEASE Genes___________#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Drug perturbed Genes (LINCS data) and Disease Genes 
#' from the above step (DEGs and GWAS genes overlap)
#' Output would be the disease specific genes which are also perturbed by the drugs.

# should make a separate script with all smal code snippets used for creating these small files 
load("./data/gene.id.entrez.RData")

# get LINCS dataset
load("./data/L1000.RData")
# harmonizome = unique(fread("./data/harmonizome.tsv"))
L1000 = L1000[,c(5,1,7)]
length(unique(L1000$chembl.id))
length(unique(L1000$ensembl.id))
L1000 = L1000[order(ensembl.id,decreasing = TRUE), ]
L1000 = L1000[!duplicated(L1000[,c('ensembl.id', 'chembl.id')]),] # remove suplicate entries
L1000 = merge(L1000,gene.id.entrez,by="ensembl.id") # merging with ENTREZ ID, since we need only ones with ENTREZ IDs for SPIA calculation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~Optional~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' In our case we want to investigate all drugs which went until at least any clinical trial phases
#' if users want to investigate all drugs they can comment out next 4 lines

load("./data/drug2715details.RData")
dmap = data.table(dmap[,c(1,2,3)])
dmap = dmap[phase==4|phase==3|phase==2|phase==1] #673 Drugs
L1000 = merge(dmap[,c(1,2)],L1000,by="chembl.id")
length(unique(L1000$chembl.id))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##### create Disease-Genes list #####
load("./data/disease.genes50.RData")
# disease.genes = fread("./data/disease.genes50.tsv")
# considering all those gene sets where DEGs and GWAS came from same diseases.
DisGen = disease.genes[same.disease == TRUE & overlap > 0] 
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.DEGs),] #remove duplicated rows based on one column
DisGenX = unique(DisGen %>% 
                         mutate(commonGenes = strsplit(as.character(commonGenes), ",")) %>% 
                         unnest(commonGenes))
DisGenX = DisGenX[,c(1,2,13)]
# cleaning up symbols
DisGenX$commonGenes = gsub("^c\\(", "", DisGenX$commonGenes)
DisGenX$commonGenes = gsub("\\)", "", DisGenX$commonGenes)
DisGenX$commonGenes = gsub("\"", "", DisGenX$commonGenes)
DisGenX$commonGenes = trimws(DisGenX$commonGenes, which = "both")

DisGen.list = split(DisGenX, DisGenX$efo.id.DEGs)

# create vector with all disease, ensembl and chemical IDs
efo.ids = unique(DisGen$efo.id.DEGs)
ensembl.ids = unique(keys(EnsDb.Hsapiens.v86))
chembl.ids = unique(L1000[, chembl.id])
load("./data/drug2disease.RData")
opentargets.drugs = unique(fread("./data/opentargets.drugs.tsv"))

## Signifcant overlap calculation
drugPdisease.genes <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
  this.efo.id = efo.ids[i]
  # loop through drugs
  foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
    this.chembl.id = chembl.ids[j]
    degs = unique(L1000[chembl.id == this.chembl.id, ensembl.id])
    gags = DisGen.list[[this.efo.id]]$commonGenes
    if (length(degs) > 0 && length(gags) > 0) {
      # test significance
      tmp = SignificantOverlap(degs, gags, ensembl.ids)
      # annotate
      tmp = cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
      tmp = merge(unique(drug2disease[, .(efo.id, efo.term)]), tmp, by = "efo.id")
      tmp = merge(unique(L1000[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
      # add dummy variable
      tmp[, existing.indication := ifelse(nrow(drug2disease[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
      setcolorder(tmp, c("chembl.id","chembl.name", "efo.id", "efo.term", "DEGs", "GWASs", "overlap", "universe","commonGenes", "odds.ratio", "p.value", "existing.indication"))
    }
  }
}

drugPdisease.genes = unique(drugPdisease.genes)
drugPdisease.genes = drugPdisease.genes[overlap > 0]
length(unique(drugPdisease.genes$chembl.id))
length(unique(drugPdisease.genes$efo.id))
save(drugPdisease.genes,file="./data/drugPdisease.genes50.RData")

#load("./data/drugPdisease.genes50.RData")
# correct p-values
drugPdisease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
drugPdisease.genes = drugPdisease.genes[p.adjusted < 0.05]
save(drugPdisease.genes,file="./data/drugPdisease.genes50.padj.RData")
drugPdisease.genes = drugPdisease.genes[p.adjusted < 1e-05]
save(drugPdisease.genes,file="./data/drugPdisease.genes50.padj1e-5.RData")
#md2 = unique(min.drugs[,c(1,3,12)]) #filtering columns for merging indication area to our super drugs
#drugPdisease.genes = merge(drugPdisease.genes,md2,by=c("chembl.id","efo.id"))
drugPdisease.genes = drugPdisease.genes[order(p.adjusted), ]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____Drug perturbed Genes to GWAS Genes_______####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.
#' 

# GWASs = unique(fread("./data/stopgap.tsv"))
load("./data/drug2disease.RData")
load("./data/GWASs.RData")
GWASs = data.table(GWASs[,c(4,5,2,3)])
tmp = dplyr::count(GWASs, efo.id)

# disgen_plot1 = ggplot(tmp, aes(x=n)) + geom_freqpoly(color="darkred",bins = 30)+ ggtitle("Genes per Disease frequency plot")
# 
tmp2 = subset(tmp, n >= 50)
# disgen_plot2 = ggplot(tmp2, aes(x=n)) + geom_freqpoly(color="darkred", bins = 30)+ ggtitle("Genes (>= 50) per Disease frequency plot")
# 
tmp3 = subset(tmp, n >= 100)
# disgen_plot3 = ggplot(tmp3, aes(x=n)) + geom_freqpoly(color="darkred", bins = 30)+ ggtitle("Genes (>= 100) per Disease frequency plot")
# 
tmp4 = subset(tmp, n >= 500)
# disgen_plot4 = ggplot(tmp4, aes(x=n)) + geom_freqpoly(color="darkred", bins = 30)+ ggtitle("Genes (>= 500) per Disease frequency plot")
# 
# jpeg(file="./data/frequency_plots_GWASs_diseaseGenes.jpeg", width=1800, height=1980, res=200)
# plot_grid(disgen_plot1,disgen_plot2, disgen_plot3,disgen_plot4,align = "h", ncol = 2,nrow =2, rel_heights = c(1/4, 1/4))
# dev.off()


GWASs = merge(tmp2,GWASs, by="efo.id")
GWASs = data.table(GWASs[,c(1,3,4,5)])

tmp = dplyr::count(L1000, chembl.id)
tmp2 = subset(tmp, n >= 10) #we are selecting drugs that alter expresiion of at least 10 genes
L1000 = merge(L1000,tmp2,by="chembl.id")
L1000 = L1000[,c(1,2,3,4,5)]
#update chembl & efo.ids again
chembl.ids = unique(L1000[, chembl.id])
efo.ids = unique(GWASs$efo.id)
GWASs.list = split(GWASs, GWASs$efo.id)
# create gene universes for Fisher's test
ensembl.ids = unique(keys(EnsDb.Hsapiens.v86))

## Signifcant overlap calculation
drugGWAS.genes <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
  this.efo.id = efo.ids[i]
  # loop through drugs
  foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
    this.chembl.id = chembl.ids[j]
    degs = unique(L1000[chembl.id == this.chembl.id, ensembl.id])
    gags = GWASs.list[[this.efo.id]]$ensembl.id
    if (length(degs) > 0 && length(gags) > 0) {
      # test significance
      tmp = SignificantOverlap(degs, gags, ensembl.ids)
      # annotate
      tmp = cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
      tmp = merge(unique(GWASs[, .(efo.id, efo.term)]), tmp, by = "efo.id")
      tmp = merge(unique(L1000[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
      # add dummy variable
      tmp[, existing.indication := ifelse(nrow(drug2disease[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
      setcolorder(tmp, c("chembl.id", "chembl.name","efo.id", "efo.term", "DEGs", "GWASs", "overlap", "universe","commonGenes", "odds.ratio", "p.value", "existing.indication"))
    }
  }
}

drugGWAS.genes = unique(drugGWAS.genes)

# correct p-values
drugGWAS.genes[p.value == 0, p.value := 3e-324]
drugGWAS.genes = drugGWAS.genes[overlap > 0]
drugGWAS.genes = drugGWAS.genes[order(p.value), ]
drugGWAS.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
save(drugGWAS.genes,file="./data/drugGWAS.genes50.RData")

drugGWAS.genes = drugGWAS.genes[p.adjusted < 1e-05]
save(drugGWAS.genes,file="./data/drugGWAS.genes50.padj1e-5.RData")

drugGWAS.genes = drugGWAS.genes[p.adjusted < 1e-10]
save(drugGWAS.genes,file="./data/drugGWAS.genes50.padj1e-10.RData")


load("./data/drugGWAS.genes50.padj1e-5.RData")
load("./data/drugGWAS.genes50.padj1e-10.RData")
length(unique(drugGWAS.genes$chembl.id))
length(unique(drugGWAS.genes$efo.id))

load("./data/drugPdisease.genes.48D.padj1e-5.RData")
length(unique(drugPdisease.genes$efo.id))
