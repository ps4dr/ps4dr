library(EnsDb.Hsapiens.v86)

library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

library(SPIA)

hgnc = as.character(unique(gene.id$ENTREZ))

lfc.df = data.table(lfc.efo[["EFO_0000249"]])
gene.ad1 = unique(lfc.df$ENTREZ)


SignificantOverlap <- function(set1, set2, universe) {
  set1 <- unique(set1)
  set2 <- unique(set2)
  x <- list(intersect(set1,set2))
  a <- sum(set1 %in% set2)
  b <- length(set1) - a
  c <- length(set2) - a
  d <- length(universe) - a - b - c
  fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  data.table(Pathway.Genes = length(set1), GWAGs = length(set2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}



fisher.kegg.ad <- foreach (i = seq(kegg.list), .combine = rbind, .errorhandling = "remove") %do% {
  kegg.mechanisms <- names(kegg.list)[i]
  tmp <- SignificantOverlap(kegg.list[[i]], gene.ad1, hgnc)
  
  count1=0
  count2=0
  for(gene1 in 1:lengths(tmp$commonGenes)){
    tm = tmp$commonGenes[[1]][[gene1]]
    for(gene2 in 1:nrow(lfc.df)) {
      tm2 = lfc.df[gene2][[5]]
      if(tm == tm2 & lfc.df[gene2][[4]] >= 0){
        count1 = count1 +1
      }
      else if(tm == tm2 & lfc.df[gene2][[4]] < 0){
        count2 = count2 +1
      }
    }
    if (count1 > count2){
      tmp[, status := 'Activated']
    }
    if (count1 < count2){
      tmp[, status := 'Inhibited']
    }
    if (count1 == count2){
      tmp[, status := 'Neutral']
    }
  }
  tmp <- cbind(kegg.mechanisms, tmp)
}

fisher.kegg.ad <- fisher.kegg.ad[order(p.value), ]
fisher.kegg.ad[, p.adjusted := p.adjust(p.value, method = "fdr")]
fisher.kegg.ad = fisher.kegg.ad[p.adjusted < 0.05]




#####_____Datasets________#####
load("./data/disease.genes50.lfc.RData")

#################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____DEGs to GWAS Genes >>> DISEASE Genes_____####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.

# create gene universes for Fisher's test
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))
load("./data/geneID.RData")
hgnc.symbols <- unique(gene.id$hgnc_symbol)

# split data table by disease
GWAGs.list <- split(GWAGs, GWAGs$efo.id)
DEGs.list <- split(DEGs, DEGs$efo.id)

## Signifcant overlap calculation
disease.genes <- foreach (i = seq(DEGs.list), .combine = rbind, .errorhandling = "remove") %do% {
  efo.id.DEGs <- names(DEGs.list)[i]
  foreach (j = seq(GWAGs.list), .combine = rbind, .errorhandling = "remove") %dopar% {
    efo.id.GWAGs <- names(GWAGs.list)[j]
    tmp <- SignificantOverlap(DEGs.list[[efo.id.DEGs]]$ensembl.id, GWAGs.list[[efo.id.GWAGs]]$ensembl.id, ensembl.ids)
    tmp <- cbind(efo.id.DEGs, efo.id.GWAGs, tmp)
    tmp <- merge(tmp, unique(GWAGs[, .(efo.id, efo.term)]), by.x = "efo.id.GWAGs", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
    tmp <- merge(tmp, unique(DEGs[, .(efo.id, efo.term)]), by.x = "efo.id.DEGs", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".GWAGs", ".DEGs"))
    setcolorder(tmp, c("efo.id.DEGs", "efo.term.DEGs", "efo.id.GWAGs", "efo.term.GWAGs", "DEGs", "GWAGs", "overlap", "universe","commonGenes", "odds.ratio", "p.value"))
  }
}

lfc.efo[[efo.id.DEGs]]$ensembl.id
lfc.efo[[i]][["ensembl.id"]]

fisher.kegg.ad <- foreach (i = seq(lfc.efo), .combine = rbind, .errorhandling = "remove") %do% {
  efo.id <- names(lfc.efo)[i]
  foreach (j = seq(kegg.list), .combine = rbind, .errorhandling = "remove") %dopar% {
    kegg.mechanisms <- names(kegg.list)[j]
    tmp <- SignificantOverlap(lfc.efo[[i]][["ENTREZ"]],kegg.list[[j]], hgnc.symbols)
    
    count1=0
    count2=0
    for(gene1 in 1:lengths(tmp$commonGenes)){
      tm = tmp$commonGenes[[1]][[gene1]]
      for(gene2 in 1:nrow(lfc.df)) {
        tm2 = lfc.df[gene2][[5]]
        if(tm == tm2 & lfc.df[gene2][[4]] >= 0){
          count1 = count1 +1
        }
        else if(tm == tm2 & lfc.df[gene2][[4]] < 0){
          count2 = count2 +1
        }
      }
      if (count1 > count2){
        tmp[, status := 'Activated']
      }
      if (count1 < count2){
        tmp[, status := 'Inhibited']
      }
      if (count1 == count2){
        tmp[, status := 'Neutral']
      }
    }
    tmp <- cbind(kegg.mechanisms, tmp)
  }
}





disease.genes <- unique(disease.genes)
# correct p-values
disease.genes[p.value == 0, p.value := 3e-324]
disease.genes <- disease.genes[order(p.value), ]
disease.genes = disease.genes[overlap > 0]
disease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
disease.genes <- fread("./data/disease.genes50.tsv")
disease.genes = disease.genes[p.adjusted <= 0.05]
# add dummy variable
disease.genes[, same.disease := ifelse(efo.id.DEGs == efo.id.GWAGs, TRUE, FALSE)]
#fwrite(disease.genes, "./data/disease.genes50.tsv", sep = "\t")

