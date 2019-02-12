#' 1st script
#' summary:
#' cleaning up STOPGAP data



library(EnsDb.Hsapiens.v86)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)
library(BiocParallel)
register(MulticoreParam(workers = 1))
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)
library(biomaRt)
setwd("/home/memon/projects/msdrp/")

# function to calculate enrichment
checkOverlapSignificance <- function(set1, set2, universe) {
  set1 <- unique(set1)
  set2 <- unique(set2)
  a <- sum(set1 %in% set2)
  x <- list(intersect(set1,set2))
  b <- length(set1) - a
  c <- length(set2) - a
  d <- length(universe) - a - b - c
  fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  data.table(Drug.perturbed.Genes = length(set1), Disease.Genes = length(set2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}


# read datasets
fisher.genes <- fread("./data/fisher.genes.commongenes.tsv")
fisher.drugs <- fread("./data/fisher.drugs.commongenes.tsv")
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))

##### create Drug-Genes list #####
# small.drugs = fisher.drugs[efo.id == 'EFO_0000249']
# small.drugs <- small.drugs[p.adjusted <= 0.05]

min.drugs <- fisher.drugs[p.adjusted <= 0.05]
length(unique(fisher.drugs$chembl.id))
length(unique(min.drugs$chembl.id))
# split column with multiple elements into multiple rows
min.drugsX <- unique(min.drugs %>% 
                        mutate(commonGenes = strsplit(as.character(commonGenes), "\\|")) %>% 
                        unnest(commonGenes))
min.drugs.list <- split(min.drugsX, list(min.drugsX$chembl.id,min.drugsX$efo.id)) #
min.drugs.list= Filter(function(x) dim(x)[1] > 0, min.drugs.list) # remove empty lists

#chembl.ids <- unique(small.drugs$chembl.id)
chembl.ids <- unique(names(min.drugs.list))

##### create Disease-Genes list #####
small.genes = fisher.genes[same.disease == TRUE & overlap > 0]
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
small.genes = small.genes[!duplicated(small.genes$efo.id.opentargets.degs),] #remove duplicated rows based on one column
small.genesX <- unique(small.genes %>% 
                          mutate(commonGenes = strsplit(as.character(commonGenes), "\\|")) %>% 
                          unnest(commonGenes))
small.genes.list <- split(small.genesX, small.genesX$efo.id.opentargets.degs)
efo.ids <- unique(small.genes$efo.id.opentargets.degs)

##### Find Drugs which influences Disease Specific Genes #####

system.time(
  super.drugs <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
    efo.id <- efo.ids[i]
    # loop through drugs
    foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
      chembl.id <- chembl.ids[j]
      if (small.genes.list[[efo.id]][["efo.id.stopgap"]][[1]] == min.drugs.list[[chembl.id]][["efo.id"]][[1]]) {
        tmp <- checkOverlapSignificance(min.drugs.list[[chembl.id]]$commonGenes, small.genes.list[[efo.id]]$commonGenes, ensembl.ids)
        chembl.id <- gsub("\\..*","",chembl.id)
        tmp <- cbind(chembl.id, efo.id, tmp)
        tmp <- merge(tmp, unique(small.genes[, .(efo.id.opentargets.degs, efo.term.opentargets.degs)]), by.x = "efo.id", by.y = "efo.id.opentargets.degs", all.x = TRUE, all.y = FALSE)
        tmp <- merge(tmp, unique(min.drugs[, .(chembl.id, chembl.name)]), by = "chembl.id", all.x = TRUE, all.y = FALSE)
        setcolorder(tmp, c("chembl.id", "chembl.name", "efo.id", efo.term="efo.term.opentargets.degs", "Drug.perturbed.Genes","Disease.Genes", "overlap", "universe","commonGenes", "odds.ratio", "p.value")) 
      }
    }
  }
)
save(super.drugs,file="./data/super.drugs.RData")

load("./data/super.drugs.RData")

# correct p-values
super.drugs[, p.adjusted := p.adjust(p.value, method = "fdr")]
super.drugs <- super.drugs[p.adjusted <= 0.05]
md2 <- unique(min.drugs[,c(1,3,12)]) #filtering columns for merging indication area to our super drugs
super.drugs <- merge(super.drugs,md2,by=c("chembl.id","efo.id"))
super.drugs <- super.drugs[order(p.adjusted), ]

#save(super.drugs,file="./data/super.drugs.RData")
load(("./data/super.drugs.RData"))
super.drugs <- super.drugs[p.adjusted <= 1e-5]


#### Get Gene Symbol for Ensembl IDs for hgnc2ensembl ID mapping #####

ensembl = useEnsembl(biomart="ensembl", version=92, dataset="hsapiens_gene_ensembl")
val <- c(1:23,"X","Y")
chr_genes <- getBM(attributes=c('chromosome_name','hgnc_symbol','ensembl_gene_id'), 
                   filters ='chromosome_name', values =val, mart = ensembl)

chr_genes$hgnc_symbol <- gsub("^$", NA, chr_genes$hgnc_symbol)
chr_genes <- chr_genes[which(!is.na(chr_genes$hgnc_symbol)),]

save(chr_genes,file="./data/hgnc2ensembl.RData")
##################################################




#kegg <- fread("./data/kegg_gene_sets.csv",strip.white=TRUE)
#replace hgnc symbol with ensembl id
for(i in 1:ncol(kegg)){
  kegg[[i]] <- chr_genes[[2]][match(kegg[[i]], chr_genes[[1]])]
}

associated genes (GWAS)associated genes (GWAS)


