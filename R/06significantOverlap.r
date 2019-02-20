#' 6th script
#' summary:
#' calculate significant overlap between/among different data sets 
#' e.g., Drug Perturbed Genes & DEGs data,
#' Drug Perturbed Genes & GWAS data,
#' GWAS & DEGs data,
#' Drug Perturbed Genes & DEGs & GWAS data

library(EnsDb.Hsapiens.v86)

library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

library(data.table)
library(dplyr)
library(tidyr)

setwd("/home/memon/projects/msdrp/")

# enrichment calculation function for two gene sets using Fisher's exact test 
SignificantOverlap <- function(set1, set2, universe) {
  set1 <- unique(set1)
  set2 <- unique(set2)
  a <- sum(set1 %in% set2)
  x <- list(intersect(set1,set2))
  b <- length(set1) - a
  c <- length(set2) - a
  d <- length(universe) - a - b - c
  fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  data.table(DEGs = length(set1), GWAGs = length(set2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}

#####_____Datasets________#####
stopgap <- unique(fread("./data/stopgap.tsv"))
opentargets.degs <- unique(fread("./data/opentargets.degs.tsv"))


#################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____DEGs to GWAS Genes >>> DISEASE Genes_____####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.

# create gene universes for Fisher's test
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))

# split data table by disease
stopgap.list <- split(stopgap, stopgap$efo.id)
opentargets.degs.list <- split(opentargets.degs, opentargets.degs$efo.id)

## Signifcant overlap calculation
disease.genes <- foreach (i = seq(opentargets.degs.list), .combine = rbind, .errorhandling = "remove") %do% {
  efo.id.opentargets.degs <- names(opentargets.degs.list)[i]
  foreach (j = seq(stopgap.list), .combine = rbind, .errorhandling = "remove") %dopar% {
    efo.id.stopgap <- names(stopgap.list)[j]
    tmp <- SignificantOverlap(opentargets.degs.list[[efo.id.opentargets.degs]]$ensembl.id, stopgap.list[[efo.id.stopgap]]$ensembl.id, ensembl.ids)
    tmp <- cbind(efo.id.opentargets.degs, efo.id.stopgap, tmp)
    tmp <- merge(tmp, unique(stopgap[, .(efo.id, efo.term)]), by.x = "efo.id.stopgap", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
    tmp <- merge(tmp, unique(opentargets.degs[, .(efo.id, efo.term)]), by.x = "efo.id.opentargets.degs", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".stopgap", ".opentargets.degs"))
    setcolorder(tmp, c("efo.id.opentargets.degs", "efo.term.opentargets.degs", "efo.id.stopgap", "efo.term.stopgap", "DEGs", "GWAGs", "overlap", "universe","commonGenes", "odds.ratio", "p.value"))
  }
}

disease.genes <- unique(disease.genes)
# correct p-values
disease.genes[p.value == 0, p.value := 3e-324]
disease.genes <- disease.genes[order(p.value), ]
disease.genes = disease.genes[overlap > 0]
disease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
# add dummy variable
disease.genes[, same.disease := ifelse(efo.id.opentargets.degs == efo.id.stopgap, TRUE, FALSE)]
fwrite(disease.genes, "./data/disease.genes.tsv", sep = "\t")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##___________Drugs to DISEASE Genes___________#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Drug perturbed Genes (LINCS data) adn Disease Genes 
#' from the above step (DEGs and GWAS genes overlap)
#' Output would be the disease specific genes which are also perturbed by the drugs.

# should make a separate script with all smal code snippets used for creating these small files 
load("./data/gene.id.entrez.RData")

# get LINCS dataset
harmonizome <- unique(fread("./data/harmonizome.tsv"))
harmonizome = harmonizome[,c(5,1,7)]
# length(unique(harmonizome$chembl.id))
# length(unique(harmonizome$ensembl.id))
harmonizome = harmonizome[order(ensembl.id,decreasing = TRUE), ]
harmonizome = harmonizome[!duplicated(harmonizome[,c('ensembl.id', 'chembl.id')]),] # remove suplicate entries
harmonizome = merge(harmonizome,gene.id.entrez,by="ensembl.id") # merging with ENTREZ ID, since we need only ones with ENTREZ IDs for SPIA calculation

#' In our case we want to investigate all drugs which went until at least any clinical trial phases
#' if users want to investigate all drugs they can comment out next 4 lines
load("./data/drug2715details.RData")
dmap = data.table(dmap[,c(1,2,3)])
dmap = dmap[phase==4|phase==3|phase==2|phase==1]
harmonizome = merge(dmap[,c(1,2)],harmonizome,by="chembl.id")

##### create Disease-Genes list #####
disease.genes <- fread("./data/disease.genes.tsv")
# considering all those gene sets where DEGs and GWAS came from same diseases.
DisGen = disease.genes[same.disease == TRUE & overlap > 0] 
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.opentargets.degs),] #remove duplicated rows based on one column
DisGenX <- unique(DisGen %>% 
                         mutate(commonGenes = strsplit(as.character(commonGenes), "\\|")) %>% 
                         unnest(commonGenes))
DisGen.list <- split(DisGenX, DisGenX$efo.id.opentargets.degs)

# create vector with all disease, ensembl and chemical IDs
efo.ids <- unique(DisGen$efo.id.opentargets.degs)
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))
chembl.ids <- unique(harmonizome[, chembl.id])
opentargets.drugs <- unique(fread("./data/opentargets.drugs.tsv"))

## Signifcant overlap calculation
drugPdisease.genes <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
  this.efo.id <- efo.ids[i]
  # loop through drugs
  foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
    this.chembl.id <- chembl.ids[j]
    degs <- unique(harmonizome[chembl.id == this.chembl.id, ensembl.id])
    gags <- DisGen.list[[this.efo.id]]$commonGenes
    if (length(degs) > 0 && length(gags) > 0) {
      # test significance
      tmp <- SignificantOverlap(degs, gags, ensembl.ids)
      # annotate
      tmp <- cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
      tmp <- merge(unique(opentargets.drugs[, .(efo.id, efo.term)]), tmp, by = "efo.id")
      tmp <- merge(unique(harmonizome[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
      # add dummy variable
      tmp[, existing.indication := ifelse(nrow(opentargets.drugs[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
      setcolorder(tmp, c("chembl.id","chembl.name", "efo.id", "efo.term", "DEGs", "GWAGs", "overlap", "universe","commonGenes", "odds.ratio", "p.value", "existing.indication"))
    }
  }
}

drugPdisease.genes <- unique(drugPdisease.genes)
drugPdisease.genes = drugPdisease.genes[overlap > 0]
length(unique(drugPdisease.genes$chembl.id))
length(unique(drugPdisease.genes$efo.id))
#save(drugPdisease.genes,file="./data/drugPdisease.genes.48D.RData")

load("./data/drugPdisease.genes.48D.RData")
# correct p-values
drugPdisease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
drugPdisease.genes <- drugPdisease.genes[p.adjusted < 0.05]
#save(drugPdisease.genes,file="./data/drugPdisease.genes.48D.padj.RData")
drugPdisease.genes <- drugPdisease.genes[p.adjusted < 1e-05]
save(drugPdisease.genes,file="./data/drugPdisease.genes.48D.padj1e-5.RData")
#md2 <- unique(min.drugs[,c(1,3,12)]) #filtering columns for merging indication area to our super drugs
#drugPdisease.genes <- merge(drugPdisease.genes,md2,by=c("chembl.id","efo.id"))
drugPdisease.genes <- drugPdisease.genes[order(p.adjusted), ]
load("./data/drug2715details.RData")
dmap = data.table(dmap[,c(1,2,3)])
dmap = dmap[phase == 4]
dmap = data.table(dmap[,c(1,2)])

drugPdisease.genes <- merge(dmap,drugPdisease.genes,by="chembl.id")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____Drug perturbed Genes to GWAS Genes_______####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.

#update chembl & efo.ids again
chembl.ids <- unique(harmonizome[, chembl.id])
efo.ids <- unique(stopgap$efo.id)

## Signifcant overlap calculation
drugGWAS.genes <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
  this.efo.id <- efo.ids[i]
  # loop through drugs
  foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
    this.chembl.id <- chembl.ids[j]
    degs <- unique(harmonizome[chembl.id == this.chembl.id, ensembl.id])
    gags <- stopgap.list[[this.efo.id]]$ensembl.id
    if (length(degs) > 0 && length(gags) > 0) {
      # test significance
      tmp <- SignificantOverlap(degs, gags, ensembl.ids)
      # annotate
      tmp <- cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
      tmp <- merge(unique(stopgap[, .(efo.id, efo.term)]), tmp, by = "efo.id")
      tmp <- merge(unique(harmonizome[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
      # add dummy variable
      tmp[, existing.indication := ifelse(nrow(opentargets.drugs[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
      setcolorder(tmp, c("chembl.id", "chembl.name","efo.id", "efo.term", "DEGs", "GWAGs", "overlap", "universe","commonGenes", "odds.ratio", "p.value", "existing.indication"))
    }
  }
}

drugGWAS.genes <- unique(drugGWAS.genes)

# correct p-values
drugGWAS.genes[p.value == 0, p.value := 3e-324]
drugGWAS.genes = drugGWAS.genes[overlap > 0]
drugGWAS.genes <- drugGWAS.genes[order(p.value), ]
drugGWAS.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
fwrite(drugGWAS.genes, "./data/drugGWAS.genes.tsv", sep = "\t")

drugGWAS.genes <- drugGWAS.genes[p.adjusted < 1e-05]
save(drugGWAS.genes,file="./data/drugGWAS.genes.48D.padj1e-5.RData")

drugGWAS.genes <- drugGWAS.genes[p.adjusted < 1e-10]
save(drugGWAS.genes,file="./data/drugGWAS.genes.48D.padj1e-10.RData")

length(unique(drugGWAS.genes$chembl.id))
length(unique(drugGWAS.genes$efo.id))
