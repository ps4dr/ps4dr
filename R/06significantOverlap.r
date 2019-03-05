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
library(ggplot2)
library(cowplot)

setwd("/home/memon/projects/msdrp/")

# enrichment calculation function for two gene sets using Fisher's exact test 
SignificantOverlap <- function(set1, set2, universe) {
  set1 <- unique(set1)
  set2 <- unique(set2)
  x <- list(intersect(set1,set2))
  a <- sum(set1 %in% set2)
  b <- length(set1) - a
  c <- length(set2) - a
  d <- length(universe) - a - b - c
  fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="greater")
  data.table(DEGs = length(set1), GWAGs = length(set2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}

#####_____Datasets________#####
GWAGs <- unique(fread("./data/stopgap.tsv"))
tmp <- dplyr::count(GWAGs, efo.id)
tmp2 <- subset(tmp, n > 50)
GWAGs = merge(GWAGs,tmp2,by="efo.id")

DEGs <- unique(fread("./data/opentargets.degs.tsv"))
tmp <- dplyr::count(DEGs, efo.id)
tmp2 <- subset(tmp, n > 50)
DEGs = merge(DEGs,tmp2,by="efo.id")

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
fwrite(disease.genes, "./data/disease.genes50.tsv", sep = "\t")

names(disease.genes)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##___________Drugs to DISEASE Genes___________#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Drug perturbed Genes (LINCS data) and Disease Genes 
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##__________________Optional__________________#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' In our case we want to investigate all drugs which went until at least any clinical trial phases
#' if users want to investigate all drugs they can comment out next 4 lines
load("./data/drug2715details.RData")
dmap = data.table(dmap[,c(1,2,3)])
dmap = dmap[phase==4|phase==3|phase==2|phase==1]
harmonizome = merge(dmap[,c(1,2)],harmonizome,by="chembl.id")

##### create Disease-Genes list #####
disease.genes <- fread("./data/disease.genes50.tsv")
# considering all those gene sets where DEGs and GWAS came from same diseases.
DisGen = disease.genes[same.disease == TRUE & overlap > 0] 
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.DEGs),] #remove duplicated rows based on one column
DisGenX <- unique(DisGen %>% 
                         mutate(commonGenes = strsplit(as.character(commonGenes), "\\|")) %>% 
                         unnest(commonGenes))
DisGen.list <- split(DisGenX, DisGenX$efo.id.DEGs)

# create vector with all disease, ensembl and chemical IDs
efo.ids <- unique(DisGen$efo.id.DEGs)
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
save(drugPdisease.genes,file="./data/drugPdisease.genes50.RData")

#load("./data/drugPdisease.genes50.RData")
# correct p-values
drugPdisease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
drugPdisease.genes <- drugPdisease.genes[p.adjusted < 0.05]
#save(drugPdisease.genes,file="./data/drugPdisease.genes.48D.padj.RData")
drugPdisease.genes <- drugPdisease.genes[p.adjusted < 1e-05]
save(drugPdisease.genes,file="./data/drugPdisease.genes50.padj1e-5.RData")
#md2 <- unique(min.drugs[,c(1,3,12)]) #filtering columns for merging indication area to our super drugs
#drugPdisease.genes <- merge(drugPdisease.genes,md2,by=c("chembl.id","efo.id"))
drugPdisease.genes <- drugPdisease.genes[order(p.adjusted), ]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____Drug perturbed Genes to GWAS Genes_______####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.
#' 
opentargets.drugs <- unique(fread("./data/opentargets.drugs.tsv"))
GWAGs <- unique(fread("./data/stopgap.tsv"))
GWAGs = data.table(GWAGs[,c(3,4,1,2)])
tmp <- dplyr::count(GWAGs, efo.id)

# disgen_plot1 <- ggplot(tmp, aes(x=n)) + geom_freqpoly(color="darkred",bins = 30)+ ggtitle("Genes per Disease frequency plot")
# 
tmp2 <- subset(tmp, n >= 50)
# disgen_plot2 <- ggplot(tmp2, aes(x=n)) + geom_freqpoly(color="darkred", bins = 30)+ ggtitle("Genes (>= 50) per Disease frequency plot")
# 
tmp3 <- subset(tmp, n >= 100)
# disgen_plot3 <- ggplot(tmp3, aes(x=n)) + geom_freqpoly(color="darkred", bins = 30)+ ggtitle("Genes (>= 100) per Disease frequency plot")
# 
tmp4 <- subset(tmp, n >= 500)
# disgen_plot4 <- ggplot(tmp4, aes(x=n)) + geom_freqpoly(color="darkred", bins = 30)+ ggtitle("Genes (>= 500) per Disease frequency plot")
# 
# jpeg(file="./data/frequency_plots_GWAGs_diseaseGenes.jpeg", width=1800, height=1980, res=200)
# plot_grid(disgen_plot1,disgen_plot2, disgen_plot3,disgen_plot4,align = "h", ncol = 2,nrow =2, rel_heights = c(1/4, 1/4))
# dev.off()


GWAGs = merge(tmp2,GWAGs, by="efo.id")
GWAGs = data.table(GWAGs[,c(1,3,4,5)])

tmp <- dplyr::count(harmonizome, chembl.id)
tmp2 <- subset(tmp, n > 50)
harmonizome = merge(harmonizome,tmp2,by="chembl.id")
harmonizome = harmonizome[,c(1,2,3,4,5)]
#update chembl & efo.ids again
chembl.ids <- unique(harmonizome[, chembl.id])
efo.ids <- unique(GWAGs$efo.id)
GWAGs.list <- split(GWAGs, GWAGs$efo.id)
# create gene universes for Fisher's test
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))

## Signifcant overlap calculation
drugGWAS.genes <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
  this.efo.id <- efo.ids[i]
  # loop through drugs
  foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
    this.chembl.id <- chembl.ids[j]
    degs <- unique(harmonizome[chembl.id == this.chembl.id, ensembl.id])
    gags <- GWAGs.list[[this.efo.id]]$ensembl.id
    if (length(degs) > 0 && length(gags) > 0) {
      # test significance
      tmp <- SignificantOverlap(degs, gags, ensembl.ids)
      # annotate
      tmp <- cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
      tmp <- merge(unique(GWAGs[, .(efo.id, efo.term)]), tmp, by = "efo.id")
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
save(drugGWAS.genes,file="./data/drugGWAS.genes50.RData")

drugGWAS.genes <- drugGWAS.genes[p.adjusted < 1e-05]
save(drugGWAS.genes,file="./data/drugGWAS.genes50.padj1e-5.RData")

drugGWAS.genes <- drugGWAS.genes[p.adjusted < 1e-10]
save(drugGWAS.genes,file="./data/drugGWAS.genes50.padj1e-10.RData")

length(unique(drugGWAS.genes$chembl.id))
length(unique(drugGWAS.genes$efo.id))
load("./data/drugGWAS.genes.padj1e-10.RData")
load("./data/drugPdisease.genes.48D.padj1e-5.RData")
length(unique(drugPdisease.genes$efo.id))
