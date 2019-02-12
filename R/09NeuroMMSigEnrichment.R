library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)
library(BiocParallel)
register(MulticoreParam(workers = 1))
library(ggplot2)
library(pROC)
library(dplyr)
library(EnsDb.Hsapiens.v86)

setwd("/home/memon/projects/msdrp/")

# function to calculate enrichment
checkOverlapSignificance <- function(Pathway.Genes, DEG.GWA.Genes, universe) {
  Pathway.Genes <- unique(Pathway.Genes)
  DEG.GWA.Genes <- unique(DEG.GWA.Genes)
  a <- sum(Pathway.Genes %in% DEG.GWA.Genes)
  b <- length(Pathway.Genes) - a
  c <- length(DEG.GWA.Genes) - a
  d <- length(universe) - a - b - c
  fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  data.table(Pathway.Genes = length(Pathway.Genes), DEG.GWA.Genes = length(DEG.GWA.Genes), overlap = a, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}

delete.emptyList  <-  function(x.list){ 
  x.list[unlist(lapply(x.list, length) != 0)] 
} 

# read datasets
stopgap <- unique(fread("./data/stopgap.tsv"))
#opentargets.drugs <- unique(fread("./data/opentargets.drugs.tsv"))
opentargets.degs <- unique(fread("./data/opentargets.degs.tsv"))
#harmonizome <- unique(fread("./data/harmonizome.tsv"))

# create universes for Fisher's test
ensembl.ids <- keys(EnsDb.Hsapiens.v86)

#####  Neurommsig ad gene sets  #####
#nmmsig.ad <- fread("./data/neurommsig_ad_full_genes.csv") ### enriched with neighboring genes
nmmsig.ad <- fread("./data/neurommsig_ad_genesets.csv")
nmmsig.ad.list <- as.list(nmmsig.ad)
nmmsig.ad.list <- lapply(nmmsig.ad.list, function(x) x[!x==""])

###  Neurommsig pd gene sets #####
#nmmsig.pd <- fread("./data/neurommsig_pd_full_genes.csv") ### enriched with neighboring genes
nmmsig.pd <- fread("./data/neurommsig_pd_gene_sets.csv")
nmmsig.pd.list <- as.list(nmmsig.pd)
nmmsig.pd.list <- lapply(nmmsig.pd.list, function(x) x[!x==""])

#####  KEGG gene sets  #####
kegg <- fread("./data/kegg_gene_sets.csv",strip.white=TRUE)
kegg.list <- as.list(kegg)
kegg.list <- lapply(kegg.list, function(x) x[!is.na(x)])
kegg.list <- lapply(kegg.list, function(x) x[!x == ""])
x = sapply(kegg.list, function(x) x[! length(x) > 1000])
kegg.list = delete.emptyList(x)

##### Reactome gene sets  ######
react <- fread("./data/reactome_gene_sets.csv",strip.white=TRUE)
react.list <- as.list(react)
react.list <- lapply(react.list, function(x) x[!x == ""])
x = sapply(react.list, function(x) x[! length(x) > 1000])
react.list = delete.emptyList(x)

##### WiKi Pathway gene sets  ######

wiki <- fread("./data/wikipathways_gene_sets.csv",strip.white=TRUE)
wiki.list <- as.list(wiki)
wiki.list <- lapply(wiki.list, function(x) x[!x == ""])
x = sapply(wiki.list, function(x) x[! length(x) > 1000])
wiki.list = delete.emptyList(x)


rm(kegg,react,nmmsig.ad,nmmsig.pd,x,wiki)

#####   Exploring Alzheimer's disease   #####

gwas.gene.ad <- unique(stopgap[efo.term == "alzheimer's disease"])
gwas.gene.ad = unique(gwas.gene.ad$gene.symbol)

deg.gene.ad <- unique(opentargets.degs[efo.term == "Alzheimer's disease"])
deg.gene.ad = unique(deg.gene.ad$gene.symbol)
gwas.deg.ad = (deg.gene.ad %in% gwas.gene.ad )

gene.ad = intersect(deg.gene.ad,gwas.gene.ad )
rm(gwas.deg.ad,deg.gene.ad,gwas.gene.ad)

#################################################
#' 1st script
#' summary:
#' cleaning up STOPGAP data

#####   Exploring parkinson disease   #####

gwas.gene.pd <- unique(stopgap[efo.term == "parkinson disease"])
gwas.gene.pd = unique(gwas.gene.pd$gene.symbol)

deg.gene.pd <- unique(opentargets.degs[efo.term == "Parkinson's disease"])
deg.gene.pd = unique(deg.gene.pd$gene.symbol)
gwas.deg.pd = (deg.gene.pd %in% gwas.gene.pd )

gene.pd = intersect(deg.gene.pd,gwas.gene.pd )
rm(gwas.deg.pd,deg.gene.pd,gwas.gene.pd)

#################################################

### are genes which both differentially expressed and genetically associated with disease X, enriched in each NeuroMMsig Mchanisms?

##### Fisher's test for NeuroMMSig AD Subgraphs #####

fisher.nmmsig.ad <- foreach (i = seq(nmmsig.ad.list), .combine = rbind, .errorhandling = "remove") %do% {
  ad.subgraph <- names(nmmsig.ad.list)[i]
  tmp <- checkOverlapSignificance(nmmsig.ad.list[[i]], gene.ad, ensembl.ids)
  tmp <- cbind(ad.subgraph, tmp)
}

# correct p-values
fisher.nmmsig.ad <- fisher.nmmsig.ad[order(p.value), ]
fisher.nmmsig.ad[, p.adjusted := p.adjust(p.value, method = "fdr")]
#fisher.nmmsig.ad = fisher.nmmsig.ad[p.adjusted < 0.05]
fisher.nmmsig.ad = fisher.nmmsig.ad[p.value < 0.05]

#################################

##### Fisher's test for NeuroMMSig PD Subgraphs #####

fisher.nmmsig.pd <- foreach (i = seq(nmmsig.pd.list), .combine = rbind, .errorhandling = "remove") %do% {
  pd.subgraph <- names(nmmsig.ad.list)[i]
  tmp <- checkOverlapSignificance(nmmsig.pd.list[[i]], gene.pd, ensembl.ids)
  tmp <- cbind(pd.subgraph, tmp)
}

# correct p-values
fisher.nmmsig.pd <- fisher.nmmsig.pd[order(p.value), ]
fisher.nmmsig.pd[, p.adjusted := p.adjust(p.value, method = "fdr")]
#fisher.nmmsig.pd = fisher.nmmsig.pd[p.adjusted < 0.05]
fisher.nmmsig.pd = fisher.nmmsig.pd[p.value < 0.05]

#################################

##### Enrichment of KEGG gene sets for significant AD genes  #####

fisher.kegg.ad <- foreach (i = seq(kegg.list), .combine = rbind, .errorhandling = "remove") %do% {
  kegg.mechanisms <- names(kegg.list)[i]
  tmp <- checkOverlapSignificance(kegg.list[[i]], gene.ad, ensembl.ids)
  tmp <- cbind(kegg.mechanisms, tmp)
}

fisher.kegg.ad <- fisher.kegg.ad[order(p.value), ]
fisher.kegg.ad[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.kegg.ad = fisher.kegg.ad[p.adjusted < 0.05]

#####################################

##### Enrichment of KEGG gene sets for significant AD genes  #####

fisher.kegg.pd <- foreach (i = seq(kegg.list), .combine = rbind, .errorhandling = "remove") %do% {
  kegg.mechanisms <- names(kegg.list)[i]
  tmp <- checkOverlapSignificance(kegg.list[[i]], gene.pd, ensembl.ids)
  tmp <- cbind(kegg.mechanisms, tmp)
}

fisher.kegg.pd <- fisher.kegg.pd[order(p.value), ]
fisher.kegg.pd[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.kegg.pd = fisher.kegg.pd[p.adjusted < 0.05]

#####################################

##### Enrichment of Reactome gene sets for significant AD genes  #####

fisher.react.ad <- foreach (i = seq(react.list), .combine = rbind, .errorhandling = "remove") %do% {
  react.mechanisms <- names(react.list)[i]
  tmp <- checkOverlapSignificance(react.list[[i]], gene.ad, ensembl.ids)
  tmp <- cbind(react.mechanisms, tmp)
}

fisher.react.ad <- fisher.react.ad[order(p.value), ]
fisher.react.ad[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.react.ad = fisher.react.ad[p.adjusted < 0.05]

#####################################

##### Enrichment of Reactome gene sets for significant PD genes  #####

fisher.react.pd <- foreach (i = seq(react.list), .combine = rbind, .errorhandling = "remove") %do% {
  react.mechanisms <- names(react.list)[i]
  tmp <- checkOverlapSignificance(react.list[[i]], gene.pd, ensembl.ids)
  tmp <- cbind(react.mechanisms, tmp)
}

fisher.react.pd <- fisher.react.pd[order(p.value), ]
fisher.react.pd[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.react.pd = fisher.react.pd[p.adjusted < 0.05]

#####################################

##### Enrichment of Reactome gene sets for significant AD genes  #####

fisher.wiki.ad <- foreach (i = seq(wiki.list), .combine = rbind, .errorhandling = "remove") %do% {
  wiki.mechanisms <- names(wiki.list)[i]
  tmp <- checkOverlapSignificance(wiki.list[[i]], gene.ad, ensembl.ids)
  tmp <- cbind(wiki.mechanisms, tmp)
}

fisher.wiki.ad <- fisher.wiki.ad[order(p.value), ]
fisher.wiki.ad[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.wiki.ad = fisher.wiki.ad[p.adjusted < 0.05]

#####################################

##### Enrichment of Reactome gene sets for significant AD genes  #####

fisher.wiki.pd <- foreach (i = seq(wiki.list), .combine = rbind, .errorhandling = "remove") %do% {
  wiki.mechanisms <- names(wiki.list)[i]
  tmp <- checkOverlapSignificance(wiki.list[[i]], gene.pd, ensembl.ids)
  tmp <- cbind(wiki.mechanisms, tmp)
}

fisher.wiki.pd <- fisher.wiki.pd[order(p.value), ]
fisher.wiki.pd[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.wiki.pd = fisher.wiki.pd[p.adjusted < 0.05]

#####################################

fwrite(fisher.nmmsig.ad,"./data/fisher.nmmsig.ad.padj.csv", sep = ",")
fwrite(fisher.nmmsig.pd,"./data/fisher.nmmsig.pd.padj.csv", sep = ",")

fwrite(fisher.kegg.ad,"./data/fisher.kegg.ad.padj.csv", sep = ",")
fwrite(fisher.kegg.pd,"./data/fisher.kegg.pd.padj.csv", sep = ",")

fwrite(fisher.react.ad,"./data/fisher.react.ad.padj.csv", sep = ",")
fwrite(fisher.react.pd,"./data/fisher.react.pd.padj.csv", sep = ",")

fwrite(fisher.wiki.ad,"./data/fisher.wiki.ad.padj.csv", sep = ",")
fwrite(fisher.wiki.pd,"./data/fisher.wiki.pd.padj.csv", sep = ",")
