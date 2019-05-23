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
  data.table(DEGs = length(set1), GWAGs = length(set2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}

# read datasets
stopgap <- unique(fread("./data/stopgap.tsv"))
opentargets.drugs <- unique(fread("./data/opentargets.drugs.tsv"))
opentargets.degs <- unique(fread("./data/opentargets.degs.tsv"))
harmonizome <- unique(fread("./data/harmonizome.tsv"))

#####   Exploring Alzheimer's disease   #####

alzgene.gwas <- unique(stopgap[efo.term == "alzheimer's disease"])
uniqalzgene.gwas = unique(alzgene.gwas$gene.symbol)

alzgene.deg <- unique(opentargets.degs[efo.term == "Alzheimer's disease"])
uniqalzgene.deg = unique(alzgene.deg$gene.symbol)
alz_gwas_deg = (uniqalzgene.gwas %in% uniqalzgene.deg)
table(alz_gwas_deg)
#'alz_gwas_deg
#'FALSE  TRUE 
#'2766   598 
#length(alz_gwas_deg[alz_gwas_deg == TRUE])
#' 1st script
#' summary:
#' cleaning up STOPGAP data

#################################################


# split by disease
stopgap.list <- split(stopgap, stopgap$efo.id)
opentargets.degs.list <- split(opentargets.degs, opentargets.degs$efo.id)

# create universes for Fisher's test
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))
efo.ids <- unique(c(stopgap[, efo.id], opentargets.degs[, efo.id], opentargets.drugs[, efo.id]))
chembl.ids <- unique(opentargets.drugs[, chembl.id])

#### Get Gene Symbol for Ensembl IDs 
# library(biomaRt)
# ensembl = useEnsembl(biomart="ensembl", version=93, dataset="hsapiens_gene_ensembl")
# val <- c(1:22,'X','Y','MT')
# chr_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name'), 
#                    filters ='chromosome_name', values =val, mart = ensembl)
# chr_genes$hgnc_symbol <- gsub("^$", NA, chr_genes$hgnc_symbol)
# chr_genes <- chr_genes[which(!is.na(chr_genes$hgnc_symbol)),]
# chr_genes <- unique(chr_genes)
# hgnc.list = chr_genes[[1]]
# 
# setdiff(ensembl.ids,hgnc.list)

### question 1: are genes differentially expressed in disease X enriched in genes genetically associated with disease X, compared to other diseases?
## perform Fisher's test
# loop through diseases

disease.genes <- foreach (i = seq(opentargets.degs.list), .combine = rbind, .errorhandling = "remove") %do% {
  efo.id.opentargets.degs <- names(opentargets.degs.list)[i]
  foreach (j = seq(stopgap.list), .combine = rbind, .errorhandling = "remove") %dopar% {
    efo.id.stopgap <- names(stopgap.list)[j]
    # test significance
    tmp <- SignificantOverlap(opentargets.degs.list[[efo.id.opentargets.degs]]$ensembl.id, stopgap.list[[efo.id.stopgap]]$ensembl.id, ensembl.ids)
    # annotate
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
disease.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
fish = disease.genes
disease.genes <- fread("/home/memon/projects/drug_target_prioritization/GCMap/dat/disease.genes.commongenes.tsv")
# add dummy variable
disease.genes[, same.disease := ifelse(efo.id.opentargets.degs == efo.id.stopgap, TRUE, FALSE)]

## test if distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
# adjusted p-values
mw.res <- wilcox.test(x = disease.genes[same.disease == FALSE, -log10(p.adjusted)], y = disease.genes[same.disease == TRUE, -log10(p.adjusted)])
print(mw.res)
print(mw.res$p.value)
png("./data/disease.genes.pvalues.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(disease.genes, aes(x = same.disease, y = -log10(p.adjusted))) +
    geom_boxplot(fill = "#0066ff", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(disease.genes[, -log10(p.adjusted)], c(0.03, 0.97))) +
    xlab("Same disease") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()
# odds ratios
mw.res <- wilcox.test(x = disease.genes[same.disease == FALSE, odds.ratio], y = disease.genes[same.disease == TRUE, odds.ratio]) 
print(mw.res)
print(mw.res$p.value)
png("./data/disease.genes.oddsratios.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(disease.genes, aes(x = same.disease, y = odds.ratio)) +
    geom_boxplot(fill = "#0066ff", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(disease.genes[, odds.ratio], c(0.06, 0.94))) +
    xlab("Same disease") +
    ylab("Odds ratio") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()

# perform ROC analysis
preds <- as.numeric(disease.genes[, -log10(p.adjusted)])
labls <- as.numeric(disease.genes[, same.disease])
#save(labls,preds,file="roc.RData")

#produce errors, need to do it in laptop
#roc.res <- roc(response = as.numeric(labls), predictor = as.numeric(preds), algorithm = 2, ci = TRUE, ci.method = "bootstrap", boot.n = 1000, parallel = TRUE, progress = "none")
load("/home/memon/projects/drug_target_prioritization/GCMap/roc.RData")
print(roc.res)
sp.ci <- ci.sp(roc.res, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
se.ci <- ci.se(roc.res, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
png("./data/disease.genes.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
par(pty = "s")
plot(roc.res, main = paste("AUC =", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
plot(se.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(sp.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(roc.res, add = TRUE, col = "#0066ff", lwd = 3)
dev.off()

## export
# master file
#fwrite(disease.genes, "./data/disease.genes.tsv", sep = "\t")
fwrite(disease.genes, "./data/disease.genes.commongenes.tsv", sep = "\t")
# table
fwrite(disease.genes[same.disease == TRUE & p.adjusted < 0.05, ], "./doc/TableS1.csv", sep = ",")


### question 2: are genes differentially expressed after treatment with drug for disease X enriched for genes genetically associated with disease X, compared to other diseases?

## perform Fisher's test
# loop through diseases
drugGWAS.genes <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
    this.efo.id <- efo.ids[i]
    # loop through drugs
    foreach (j = seq(chembl.ids), .combine = rbind, .errorhandling = "remove") %do% {
             this.chembl.id <- chembl.ids[j]
             degs <- unique(harmonizome[chembl.id == this.chembl.id, ensembl.id])
             gags <- stopgap.list[[this.efo.id]]$ensembl.id
             if (length(degs) > 0 && length(gags) > 0) {
                # test significance
                tmp <- checkOverlapSignificance(degs, gags, ensembl.ids)
                # annotate
                tmp <- cbind(chembl.id = this.chembl.id, efo.id = this.efo.id, tmp)
                tmp <- merge(unique(opentargets.drugs[, .(efo.id, efo.term)]), tmp, by = "efo.id")
                tmp <- merge(unique(opentargets.drugs[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
                # add dummy variable
                tmp[, existing.indication := ifelse(nrow(opentargets.drugs[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
                setcolorder(tmp, c("chembl.id", "chembl.name", "efo.id", "efo.term", "DEGs", "GWAGs", "overlap", "universe","commonGenes", "odds.ratio", "p.value", "existing.indication"))
             }
    }
}

drugGWAS.genes <- unique(drugGWAS.genes)

# correct p-values
drugGWAS.genes[p.value == 0, p.value := 3e-324]
drugGWAS.genes <- drugGWAS.genes[order(p.value), ]
drugGWAS.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]

## test if distributions significantly different with Mann-Whitney test (Wilcoxon rank sum test)
# adjusted p-values
mw.res <- wilcox.test(x = drugGWAS.genes[existing.indication == FALSE, -log10(p.adjusted)], y = drugGWAS.genes[existing.indication == TRUE, -log10(p.adjusted)])
print(mw.res)
print(mw.res$p.value)
png("./data/drugGWAS.genes.pvalues.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(drugGWAS.genes, aes(x = existing.indication, y = -log10(p.adjusted))) +
    geom_boxplot(fill = "#ff6600", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(drugGWAS.genes[, -log10(p.adjusted)], c(0.06, 0.94))) +
    xlab("Existing indication") +
    ylab("-log10(adjusted p-value)") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()
# odds ratios
mw.res <- wilcox.test(x = drugGWAS.genes[existing.indication == FALSE, odds.ratio], y = drugGWAS.genes[existing.indication == TRUE, odds.ratio]) 
print(mw.res)
print(mw.res$p.value)
png("./data/drugGWAS.genes.oddsratios.boxplots.png", width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(drugGWAS.genes, aes(x = existing.indication, y = odds.ratio)) +
    geom_boxplot(fill = "#ff6600", outlier.shape = NA) +
    coord_cartesian(ylim = quantile(drugGWAS.genes[, odds.ratio], c(0.04, 0.96))) +
    xlab("Existing indication") +
    ylab("Odds ratio") +
    theme_bw(18) +
    scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    ggtitle(paste("p-value =", sprintf("%.2e", mw.res$p.value))))
dev.off()

# perform ROC analysis
preds <- drugGWAS.genes[, -log10(p.adjusted)]
labls <- as.numeric(drugGWAS.genes[, existing.indication])
roc.res <- roc(response = labls, predictor = preds, algorithm = 2, ci = TRUE, ci.method = "bootstrap", boot.n = 1000, parallel = TRUE, progress = "none")
print(roc.res)
sp.ci <- ci.sp(roc.res, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
se.ci <- ci.se(roc.res, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
png("./data/drugGWAS.genes.roc.png", width = 6 * 150, height = 6 * 150, res = 150)
par(pty = "s")
plot(roc.res, main = paste("AUC =", round(roc.res$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
plot(se.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(sp.ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
plot(roc.res, add = TRUE, col = "#ff6600", lwd = 3)
dev.off()

## export
# master file
#fwrite(drugGWAS.genes, "./data/drugGWAS.genes.tsv", sep = "\t")
fwrite(drugGWAS.genes, "./data/drugGWAS.genes.tsv", sep = "\t")

# tables
fwrite(drugGWAS.genes[existing.indication == TRUE & p.adjusted < 0.05, ], "./doc/TableS2.csv", sep = ",")
fwrite(drugGWAS.genes[existing.indication == FALSE & p.adjusted < 1e-10, ], "./doc/TableS3.csv", sep = ",")



