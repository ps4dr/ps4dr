#' 1st script
#' summary:
#' This script calculates significant overlap between/among different data sets using Fisher's exact test
#' Three different overlaps are calculated in this script (Figure 1. "Gene set intersection")
#' 01: GWAS & DEGs data :> disease_genes,
#' 02: Drug Perturbed Genes & DEGs & GWAS data :> drugPdisease_genes,
#' 03: Drug Perturbed Genes & GWAS data :> drugGWAS_genes.

suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doParallel)))
registerDoParallel(parallel::detectCores() - 1)

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(pROC)))
suppressWarnings(suppressMessages(library(biomaRt)))

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~Enrichment Calculation Function for Two Gene Sets using Fisher's Exact Test~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SignificantOverlap <- function(geneSet1, geneSet2, universe) {
    geneSet1 = unique(geneSet1)
    geneSet2 = unique(geneSet2)
    x = list(intersect(geneSet1, geneSet2))
    a = sum(geneSet1 %in% geneSet2)
    b = length(geneSet1) - a
    c = length(geneSet2) - a
    d = length(universe) - a - b - c
    fisher = fisher.test(matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE), alternative = "greater")
    data.table(DEGs = length(geneSet1), GWASs = length(geneSet2), overlap = a, commonGenes = x, universe = length(universe), odds.ratio = fisher$estimate, p.value = fisher$p.value)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~Datasets~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"GWASs.RData"))
tmp = dplyr::count(GWASs, efo.id)
tmp2 = subset(tmp, n > 50)
GWASs = merge(GWASs, tmp2, by = "efo.id")

load(file.path(dataFolder,"DEGs.RData"))
tmp = dplyr::count(DEGs, efo.id)
tmp2 = subset(tmp, n > 50)
DEGs = merge(DEGs, tmp2, by = "efo.id")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____DEGs to GWAS Genes >>> DISEASE Genes_____####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases
#' and also recorded from GWAS.

# create gene universes for Fisher's test
load(file.path(dataFolder,"geneID_v97.RData"))
ensembl_ids = unique(gene_id$ensembl.id)

# split data table by disease
GWASs_list = split(GWASs, GWASs$efo.id)
DEGs_list = split(DEGs, DEGs$efo.id)

## Overlap Signifcance calculation
cat(sprintf("GWAS to DEGs Gene Set Overlap Signifcance Calculation\n"))
pb <- txtProgressBar(min=0, max=length(DEGs_list), style=3)

disease_genes <- foreach (i = seq(DEGs_list), .combine = rbind, .errorhandling = "remove") %do% {
    Sys.sleep(1)
    setTxtProgressBar(pb, i)
    efo.id.DEGs = names(DEGs_list)[i]
    # cat(sprintf("'GWAS 2 DEGs' overlap significance calculation for : %s, index #%d of #%d\n", efo.id.DEGs, i ,length(DEGs_list)))
    foreach (j = seq(GWASs_list), .combine = rbind, .errorhandling = "remove") %dopar% {
        efo.id.GWASs = names(GWASs_list)[j]
        tmp = SignificantOverlap(DEGs_list[[efo.id.DEGs]]$ensembl.id, GWASs_list[[efo.id.GWASs]]$ensembl.id, ensembl_ids)
        tmp = cbind(efo.id.DEGs, efo.id.GWASs, tmp)
        tmp = merge(tmp, unique(GWASs[, .(efo.id, efo.term)]), by.x = "efo.id.GWASs", by.y = "efo.id", all.x = TRUE, all.y = FALSE)
        tmp = merge(tmp, unique(DEGs[, .(efo.id, efo.term)]), by.x = "efo.id.DEGs", by.y = "efo.id", all.x = TRUE, all.y = FALSE, suffixes = c(".GWASs", ".DEGs"))
        setcolorder(tmp, c("efo.id.DEGs", "efo.term.DEGs", "efo.id.GWASs", "efo.term.GWASs", "DEGs", "GWASs", "overlap", "universe", "commonGenes", "odds.ratio", "p.value"))
    }
}
close(pb)
cat(sprintf("\n"))
# correct p-values
disease_genes[p.value == 0, p.value := 3e-324]
disease_genes = disease_genes[order(p.value),]
disease_genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
# add dummy variable
disease_genes[, same.disease := ifelse(efo.id.DEGs == efo.id.GWASs, TRUE, FALSE)]
# save(disease_genes, file = file.path(dataFolder,"results/disease_genes.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##___________Drugs to DISEASE Genes___________#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Drug perturbed Genes (LINCS data) and Disease Genes 
#' from the above step (DEGs and GWAS genes overlap)
#' Output would be the disease specific genes which are also perturbed by the drugs.

#__________map ensembl IDs to ENTREZ ID___________#
#-----get mapping among ENTREZ_HGNC_ENSEMBL_IDs---------#
load(file.path(dataFolder,"geneID_v97.RData"))
gene_id$ENTREZ = gsub("^$", NA, gene_id$ENTREZ)
gene_id = gene_id[which(! is.na(gene_id$ENTREZ)),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# get LINCS dataset
load(file.path(dataFolder,"L1000.RData"))
# harmonizome = unique(fread(file.path(dataFolder,"harmonizome.tsv")))
L1000 = L1000[, c(5, 1, 6)]
L1000 = L1000[order(ensembl.id, decreasing = TRUE),]
L1000 = L1000[!duplicated(L1000[, c('ensembl.id', 'chembl.id')]),] # remove duplicate entries
L1000 = merge(L1000, gene_id, by = "ensembl.id") # merging with ENTREZ ID, since we need only ones with ENTREZ IDs for SPIA calculation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~Optional~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' In our case we want to investigate all drugs which went until at least any clinical trial phases
#' if users want to investigate all drugs they can comment out next 4 lines

load(file.path(dataFolder,"drug2715details.RData"))
dmap = data.table(dmap[, c(1, 2, 3)])
dmap = dmap[phase == 4 | phase == 3 | phase == 2 | phase == 1] #673 Drugs
L1000 = merge(dmap[, c(1, 2)], L1000, by = "chembl.id")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##### create Disease-Genes list #####
load(file.path(dataFolder,"results/disease_genes.RData"))

# considering all those gene sets where DEGs and GWAS came from same diseases.
DisGen = disease_genes[same.disease == TRUE & overlap > 0]
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.DEGs),] #remove duplicated rows based on one column
DisGenX = unique(DisGen %>%
    mutate(commonGenes = strsplit(as.character(commonGenes), ",")) %>%
    unnest(commonGenes))
DisGenX = DisGenX[, c(1, 2, 13)]
# cleaning up symbols
DisGenX$commonGenes = gsub("^c\\(", "", DisGenX$commonGenes)
DisGenX$commonGenes = gsub("\\)", "", DisGenX$commonGenes)
DisGenX$commonGenes = gsub("\"", "", DisGenX$commonGenes)
DisGenX$commonGenes = trimws(DisGenX$commonGenes, which = "both")

DisGen.list = split(DisGenX, DisGenX$efo.id.DEGs)

# create vector with all disease, ensembl and chemical IDs
efo_ids = unique(DisGen$efo.id.DEGs)
load(file.path(dataFolder,"geneID_v97.RData"))
ensembl_ids = unique(gene_id$ensembl.id)
chembl_ids = unique(L1000[, chembl.id])
load(file.path(dataFolder,"drug2disease.RData"))

## Overlap Signifcance Calculation
cat(sprintf("Drug to Disease_Gene Set Overlap Signifcance Calculation\n"))
pb <- txtProgressBar(min=0, max=length(efo_ids), style=3)

drugPdisease_genes <- foreach (i = seq(efo_ids), .combine = rbind, .errorhandling = "remove") %do% {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  current_efo = efo_ids[i]
  # cat(sprintf("'Drug2Disease gene sets' overlap significance calculation for : %s, index #%d of #%d\n", current_efo, i ,length(efo_ids)))
  # loop through drugs
  foreach (j = seq(chembl_ids), .combine = rbind, .errorhandling = "remove") %dopar% {
      current_chembl = chembl_ids[j]
      degs = unique(L1000[chembl.id == current_chembl, ensembl.id])
      gags = DisGen.list[[current_efo]]$commonGenes
      if (length(degs) > 0 && length(gags) > 0) {
          # test significance
          tmp = as.data.table(SignificantOverlap(degs, gags, ensembl_ids))
          # annotate
          tmp = cbind(chembl.id = current_chembl, efo.id = current_efo, tmp)
          tmp = merge(unique(drug2disease[, .(efo.id, efo.term)]), tmp, by = "efo.id")
          tmp = merge(unique(L1000[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
          # add dummy variable
          tmp[, existing.indication := ifelse(nrow(drug2disease[efo.id == current_efo & chembl.id == current_chembl]) > 0, TRUE, FALSE)]
          setcolorder(tmp, c("chembl.id", "chembl.name", "efo.id", "efo.term", "DEGs", "GWASs", "overlap", "universe", "commonGenes", "odds.ratio", "p.value", "existing.indication"))
      }
  }
}
close(pb)
cat(sprintf("\n"))

# save(drugPdisease_genes, file = file.path(dataFolder,"results/drugPdisease_genes.RData"))
# 
# #load(file.path(dataFolder,"results/drugPdisease_genes50.RData"))
# # correct p-values
# drugPdisease_genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
# drugPdisease_genes = drugPdisease_genes[p.adjusted < 0.05]
# save(drugPdisease_genes, file = file.path(dataFolder,"results/drugPdisease_genes.padj.RData"))
# drugPdisease_genes = drugPdisease_genes[p.adjusted < 1e-05]
# save(drugPdisease_genes, file = file.path(dataFolder,"results/drugPdisease_genes.padj1e-5.RData"))
# #md2 = unique(min.drugs[,c(1,3,12)]) #filtering columns for merging indication area to our super drugs
# #drugPdisease_genes = merge(drugPdisease_genes,md2,by=c("chembl.id","efo.id"))


#______Following Part is optional filtering_______#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##____Drug perturbed Genes to GWAS Genes_______####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Calculating significant overlap between Differentially expressed Genes (DEGs) and GWAS genes
#' Output would be the disease specific genes which are diffrentially expressed in the diseases 
#' and also recorded from GWAS.
#' 

load(file.path(dataFolder,"drug2disease.RData"))
load(file.path(dataFolder,"GWASs.RData"))
GWASs = data.table(GWASs[, c(4, 5, 2, 3)])
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
# jpeg(file=file.path(dataFolder,"results/figures/frequency_plots_GWASs_diseaseGenes.jpeg", width=1800, height=1980, res=200))
# plot_grid(disgen_plot1,disgen_plot2, disgen_plot3,disgen_plot4,align = "h", ncol = 2,nrow =2, rel_heights = c(1/4, 1/4))
# dev.off()


GWASs = merge(tmp2, GWASs, by = "efo.id")
GWASs = data.table(GWASs[, c(1, 3, 4, 5)])

tmp = dplyr::count(L1000, chembl.id)
tmp2 = subset(tmp, n >= 10) #we are selecting drugs that alter expresiion of at least 10 genes
L1000 = merge(L1000, tmp2, by = "chembl.id")
L1000 = L1000[, c(1, 2, 3, 4, 5)]
#update chembl & efo_ids again
chembl_ids = unique(L1000[, chembl.id])
efo_ids = unique(GWASs$efo.id)
GWASs_list = split(GWASs, GWASs$efo.id)
# create gene universes for Fisher's test
load("./data/geneID_v97.RData")
ensembl_ids = unique(gene_id$ensembl.id)

## Signifcant overlap calculation
cat(sprintf("GWAS to Drug perturbed Gene Set Overlap Signifcance Calculation\n"))
pb <- txtProgressBar(min=0, max=length(efo_ids), style=3)
drugGWAS_genes <- foreach (i = seq(efo_ids), .combine = rbind, .errorhandling = "remove",.verbose = T) %dopar% {
    Sys.sleep(1)
    setTxtProgressBar(pb, i)
    current_efo = efo_ids[i]
    # loop through drugs
    foreach (j = seq(chembl_ids), .combine = rbind, .errorhandling = "remove",.verbose = T) %do% {
        current_chembl = chembl_ids[j]
        degs = unique(L1000[chembl.id == current_chembl, ensembl.id])
        gags = GWASs_list[[current_efo]]$ensembl.id
        if (length(degs) > 0 && length(gags) > 0) {
            # test significance
            tmp = SignificantOverlap(degs, gags, ensembl_ids)
            # annotate
            tmp = cbind(chembl.id = current_chembl, efo.id = current_efo, tmp)
            tmp = merge(unique(GWASs[, .(efo.id, efo.term)]), tmp, by = "efo.id")
            tmp = merge(unique(L1000[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
            # add dummy variable
            tmp[, existing.indication := ifelse(nrow(drug2disease[efo.id == current_efo & chembl.id == current_chembl]) > 0, TRUE, FALSE)]
            setcolorder(tmp, c("chembl.id", "chembl.name", "efo.id", "efo.term", "DEGs", "GWASs", "overlap", "universe", "commonGenes", "odds.ratio", "p.value", "existing.indication"))
        }
    }
}
close(pb)

# correct p-values
drugGWAS_genes[p.value == 0, p.value := 3e-324]
drugGWAS_genes = drugGWAS_genes[overlap > 0]
drugGWAS_genes = drugGWAS_genes[order(p.value),]
drugGWAS_genes[, p.adjusted := p.adjust(p.value, method = "fdr")]
save(drugGWAS_genes, file = file.path(dataFolder,"results/drugGWAS_genes.RData"))

drugGWAS_genes = drugGWAS_genes[p.adjusted < 1e-05]
save(drugGWAS_genes, file = file.path(dataFolder,"results/drugGWAS_genes.padj1e-5.RData"))

drugGWAS_genes = drugGWAS_genes[p.adjusted < 1e-10]
save(drugGWAS_genes, file = file.path(dataFolder,"results/drugGWAS_genes.padj1e-10.RData"))


# load(file.path(dataFolder,"results/drugGWAS_genes50.padj1e-5.RData"))
# load(file.path(dataFolder,"results/drugGWAS_genes50.padj1e-10.RData"))
# length(unique(drugGWAS_genes$chembl.id))
# length(unique(drugGWAS_genes$efo.id))
# 
# load(file.path(dataFolder,"results/drugPdisease_genes.48D.padj1e-5.RData"))
# length(unique(drugPdisease_genes$efo.id))
# cat(sprintf("\n"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~p-value comparison between same diseases and different diseases gene sets overlaps~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cat(sprintf("p-value comparison between same diseases and different diseases gene sets overlaps\n"))
man_wtny <- wilcox.test(x = disease_genes[same.disease == FALSE, -log10(p.adjusted)], y = disease_genes[same.disease == TRUE, -log10(p.adjusted)])
print(man_wtny)
print(man_wtny$p.value)
jpeg(file = file.path(dataFolder,"results/figures/disease_genes.pvalues.boxplots.jpeg"), width = 6 * 200, height = 6 * 150, res = 150)
print(ggplot(disease_genes, aes(x = same.disease, y = -log10(p.adjusted),fill=same.disease)) +
        geom_boxplot() +
        coord_cartesian(ylim = quantile(disease_genes[, -log10(p.adjusted)], c(0.03, 0.97))) +
        xlab("Same disease") +
        ylab("-log10(adjusted p-value") +
        theme_bw(18) +
        scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
        ggtitle(paste("p-value =", sprintf("%.2e", man_wtny$p.value))))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~ROC Curve~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# pred <- as.numeric(disease_genes[, - log10(p.adjusted)])
# resp <- as.numeric(disease_genes[, same.disease])
# 
# roc.curve <- roc(response = as.numeric(resp), predictor = as.numeric(pred), algorithm = 2, ci = TRUE, ci.method = "bootstrap", smooth = TRUE, boot.n = 1000, parallel = TRUE, progress = "none")
# print(roc.curve)
# sp_ci <- ci.sp(roc.curve, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
# se_ci <- ci.se(roc.curve, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
# jpeg(file = file.path(dataFolder,"results/figures/disease_genes.roc.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
# par(pty = "s")
# plot(roc.curve, main = paste("AUC =", round(roc.curve$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
# plot(se_ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
# plot(sp_ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
# plot(roc.curve, add = TRUE, col = "#0066ff", lwd = 3)
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
