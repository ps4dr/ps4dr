#' 2nd script
#' summary:
#' This script calculates pathway enrichment analysis using Signaling Pathway Impact (SPIA) for each of the Disease Genes 
#' that are differentially expressed and having SNPs from GWASs.
#' 01: prepare Gene Expression Data sets
#' 02: calculate SPIA

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(pROC)))
suppressWarnings(suppressMessages(library(SPIA)))
suppressWarnings(suppressMessages(library(graphite)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))

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
#~~~~~~~Prepare DEGs data set for SPIA~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"DEGs.RData"))
DEGs = DEGs[, c(1, 3, 6, 7)]

# process positive lfc
lfc_pos = DEGs[lfc >= 0]
#lfc_pos = lfc_pos[order(lfc,decreasing = TRUE), ] #lfc or pvalue we should sort the dataframe?
lfc_pos$pval.pos = abs(log10(lfc_pos$pval)) #created dummy variable for pvalue to order it decreasingly, since there are multiple lfc with same pvalue sometimes.
#lfc_pos = lfc_pos[order(pval), ]
lfc_pos = lfc_pos[order(lfc, pval.pos, decreasing = TRUE),]
lfc_pos = lfc_pos[! duplicated(lfc_pos[, c('efo.id', 'ensembl.id')]),]
#lfc_pos.efo = split(lfc_pos, lfc_pos$efo.id)
lfc_pos$pval.pos = NULL

# process negative lfc
lfc_neg = DEGs[lfc < 0]
#lfc_neg = lfc_neg[order(lfc), ] #lfc or pvalue we should sort the dataframe?
lfc_neg = lfc_neg[order(pval, lfc),]
lfc_neg = lfc_neg[! duplicated(lfc_neg[, c('efo.id', 'ensembl.id')]),]
#lfc_neg.efo = split(lfc_neg, lfc_neg$efo.id)

# combine both lfc
lfc_comb = rbind(lfc_pos, lfc_neg) #combine positive and negative lfc change back to a single data frame
lfc_comb = lfc_comb[order(pval),]
lfc_comb = lfc_comb[! duplicated(lfc_comb[, c('efo.id', 'ensembl.id')]),]
lfc_comb$pval = NULL
rm(lfc_neg, lfc_pos, DEGs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"geneID.RData"))
gene_id$entrezgene = gsub("^$", NA, gene_id$entrezgene)
gene_id = gene_id[which(! is.na(gene_id$entrezgene)),]
gene_id = data.table(gene_id)
gene_id = unique(gene_id[, c('entrezgene', 'ensembl_gene_id', 'hgnc_symbol')])
names(gene_id) = c("ENTREZ", "ensembl.id", "HGNC")
gene_id = gene_id[! duplicated(gene_id$ensembl.id),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lfc_comb = merge(lfc_comb, gene_id, by = "ensembl.id")
lfc_comb = lfc_comb[! duplicated(lfc_comb[, c('efo.id', 'ENTREZ')]),]

# lfc_efo = split(lfc_comb, lfc_comb$efo.id)
load(file.path(dataFolder,"disease_genes50.RData"))
disease_genes = disease_genes[p.adjusted <= 0.05]
DisGen = disease_genes[same.disease == TRUE &
    overlap > 0 &
    p.adjusted < 0.05]
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[! duplicated(DisGen$efo.id.DEGs),] #remove duplicated rows based on one column
DisGen = DisGen[, c(1, 2, 9)]
DisGen = unique(DisGen %>%
    mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>%
    unnest(ensembl.id))
DisGen$ensembl.id = gsub("\"", "", DisGen$ensembl.id)
DisGen$ensembl.id = gsub("c\\(", "", DisGen$ensembl.id)
DisGen$ensembl.id = gsub("\\)", "", DisGen$ensembl.id)
DisGen$ensembl.id = trimws(DisGen$ensembl.id)
DisGen$commonGenes = NULL
names(DisGen) = c("efo.id", "efo.term", "ensembl.id")
DisGen = merge(DisGen, lfc_comb, by = c('ensembl.id', 'efo.id')) # merge with harmonizome

lfc_efo = split(DisGen, DisGen$efo.term)
lfc_efo = Filter(function(x) dim(x)[1] > 10, lfc_efo) # remove diseases with very few (less than 10) genes to test

save(lfc_efo, file = file.path(dataFolder,"disease_genes50.lfc.RData"))
rm(lfc_comb)

#~~~~~~~~~~Create Named Vector for log fold changes in each disease~~~~~~~~~~#

load(file.path(dataFolder,"disease_genes50.lfc.RData"))###--get log fold change for all genes for each diseases-##

lfc_ensembl = list()
for (i in 1 : length(lfc_efo)) {
    lfc_ensembl[[i]] = setNames(as.numeric(lfc_efo[[i]][[4]]), as.character(lfc_efo[[i]][[1]]))
    names(lfc_ensembl)[[i]] = names(lfc_efo)[[i]]
}

lfc_entrez = list()
for (i in 1 : length(lfc_efo)) {
    lfc_entrez[[i]] = setNames(as.numeric(lfc_efo[[i]][[4]]), as.character(lfc_efo[[i]][[5]]))
    names(lfc_entrez)[[i]] = names(lfc_efo)[[i]]
}

lfc_entrezID = list()
for (i in 1 : length(lfc_efo)) {
    lfc_entrezID[[i]] = setNames(as.numeric(lfc_efo[[i]][[4]]), as.character(gsub("^", "ENTREZID:", lfc_efo[[i]][[5]])))
    names(lfc_entrezID)[[i]] = names(lfc_efo)[[i]]
}

lfc_hgnc = list()
for (i in 1 : length(lfc_efo)) {
    lfc_hgnc[[i]] = setNames(as.numeric(lfc_efo[[i]][[4]]), as.character(lfc_efo[[i]][[6]]))
    names(lfc_hgnc)[[i]] = names(lfc_efo)[[i]]
}

#~~~~~~Create a vector with all Gene universe to proxy Array Genes~~~~~~~#
hgnc_all = unique(gene_id$HGNC)
ensembl_all = unique(gene_id$ensembl.id)
entrez_all = unique(gene_id$ENTREZ)
entrezID_all = unique(gsub("^", "ENTREZID:", gene_id$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_hgnc, lfc_ensembl, lfc_entrez, lfc_entrezID, hgnc_all, ensembl_all, entrez_all, entrezID_all, file = file.path(dataFolder,"disease_genes50.lfc.namedVec.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~Preapre SPIA Pathway Data sets~~~~~~~~~~~#

#~~~Prepare KEGG SPIA pathway file from KEGG pathway kgml files~~~#
#' Here we assume that the user obtained the KGML (xml) files from
#' KEGG ftp site (or downloaded them one by one from the KEGG website).
#' We will not keep KEGG pathway data since the access to the KEGG ftp server requires a license

makeSPIAdata(kgml.path=file.path(dataFolder,"spia_input/kgml/"),organism="hsa",out.path=file.path(dataFolder,"spia_input/real_kegg/"))

#~~~Prepare Reactome SPIA pathway file with graphite package~~~#

# reactome <- pathways("hsapiens", "reactome")
# reactome <- convertIdentifiers(reactome, "ENTREZID")
# prepareSPIA(reactome, file.path(dataFolder,"spia_input/real_react/hsa"))

#~~~Prepare Reactome SPIA pathway file with graphite package~~~#

# biocarta <- pathways("hsapiens", "biocarta")
# biocarta <- convertIdentifiers(biocarta, "ENTREZID")
# prepareSPIA(biocarta, file.path(dataFolder,"spia_input/real_biocarta/hsa"))

##_______________________________###


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~ Disease SPIA~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' SPIA can be calcualted for KEGG, Reactome and WikiPathways
#' Here, we showed only with KEGG Pathways
#' Necessary data sets for calculating SPIA for Reactome
#' and WikiPathways are uploaded in the git repository


#~~~~~~~~~~~~~KEGG SPIA~~~~~~~~~~~~~#

load(file.path(dataFolder,"disease_genes50.lfc.namedVec.RData"))

pb <- txtProgressBar(min=0, max=length(lfc_entrez), style=3)
cat(sprintf("\n~~~~~SPIA calculation for Real KEGG Pathways~~~~~\n"))

spia_kegg = list()
for (i in 1 : length(lfc_entrez)) {
  Sys.sleep(1)
  # cat(sprintf("\n~~~~~SPIA for Disease #%d~~~~~\n", (length(lfc_entrez) + 1) -i))
  spia_kegg[[i]] = invisible(capture.output(spia(de = lfc_entrez[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")))
  setTxtProgressBar(pb, i)
}
close(pb)
names(spia_kegg) = names(lfc_entrez)

save(spia_kegg, file = file.path(dataFolder, "spia_output/spia_kegg_diseaseGenes.RData"))

# spia_kegg = lapply(spia_kegg, function(x) x[x$pNDE <= 0.05,])
# plotP(spia_kegg[[4]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~Pseudo KEGG SPIA~~~~~~~~~~~#

#' we create pseudo SPIA data sets for all three pathways
#' with random gene sets to see whether our results are by
#' random chance or meaningful indeed.


###___Create Pseudo (Random) Pathways____####

#___create universe with ENTREZID:____#
load(file.path(dataFolder,"geneID.RData"))
universe = unique(na.omit(gene_id$entrezgene)) # for kegg
# universe = unique(na.omit(gsub("^","ENTREZID:",gene_id$entrezgene))) # for reactome, biocarta

load(file.path(dataFolder,"spia_input/real_kegg/hsaSPIA.RData"))
# load(file.path(dataFolder,"spia_input/real_react/hsaSPIA.RData"))
# load(file.path(dataFolder,"spia_input/real_biocarta/hsaSPIA.RData"))

pseudo_path = path.info
rm(path.info)

#_____Shuffle first 25 matrices in pseudo_path______#
for (mylist in 1:length(pseudo_path)) {
  for (mymat in 1:(length(pseudo_path[[mylist]])-3)) { #not looping for last 3 indices since they are not matrices
    pseudo_path[[mylist]][[mymat]][] = pseudo_path[[mylist]][[mymat]][sample(length(pseudo_path[[mylist]][[mymat]]))]
  }
}

#___Replace ENTREZ ID with ENTREZ ID from Universe_____#
for (mylist in 1:length(pseudo_path)) {
  for (mymat in length(pseudo_path[[mylist]])-2) { #need to change index number to '1' for reactome and biocarta, and to '2' for kegg
    pseudo_path[[mylist]][[mymat]] = as.character(sample(universe,size = length(pseudo_path[[mylist]][[mymat]]),replace = TRUE))
  }
}

#_____change row & column names of (25) matrices______#
for (mylist in 1:length(pseudo_path)) {
  for (mymat in 1:(length(pseudo_path[[mylist]])-3)) { #not looping for last 3 indices since they are not matrices
    colnames(pseudo_path[[mylist]][[mymat]]) = pseudo_path[[mylist]][[26]] #change index to '27' for reactome and biocarta, and to '26' for kegg
    rownames(pseudo_path[[mylist]][[mymat]]) = pseudo_path[[mylist]][[26]] #change index to '27' for reactome and biocarta, and to '26' for kegg
  }
}

path.info <- pseudo_path

save(path.info,file = file.path(dataFolder,"spia_input/pseudo_kegg/hsaSPIA.RData"))
# save(path.info,file = file.path(dataFolder,"spia_input/pseudo_react/hsaSPIA.RData"))
# save(path.info,file = file.path(dataFolder,"spia_input/pseudo_biocarta/hsaSPIA.RData"))
rm(path.info,pseudo_path)
#_______________________________________________________##

#~~~~~~~~~~~~~Pseudo KEGG SPIA~~~~~~~~~~~~~#
pb <- txtProgressBar(min=0, max=length(lfc_entrez), style=3)
cat(sprintf("\n~~~~~SPIA calculation for simulated KEGG Pathways~~~~~\n"))

spia_kegg_pseudo = list()
for (i in 1 : length(lfc_entrez)) {
  Sys.sleep(1)
  # cat(sprintf("\n~~~~~SPIA (pseudo) for Disease #%d~~~~~\n", (length(lfc_entrez) + 1) -i))
  spia_kegg_pseudo[[i]] = invisible(capture.output(spia(de = lfc_entrez[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/pseudo_kegg/"), organism = "hsa")))
  setTxtProgressBar(pb, i)
}
close(pb)
names(spia_kegg_pseudo) = names(lfc_entrez)

save(spia_kegg_pseudo, file = file.path(dataFolder, "spia_output/spia_kegg_diseaseGenes_pseudo.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' delte later
# dmap = read.csv("/home/memon/projects/msdrp/2725_drugs_details.csv", header = F)
# dmap$V1 = gsub("[A-z]+\\':", "", dmap$V1)
# dmap$V1 = gsub("u", "", dmap$V1)
# dmap$V1 = gsub("\'", "", dmap$V1)
# dmap$V1 = gsub(" ", "", dmap$V1)
# dmap$V1 = gsub("\\{", "", dmap$V1)
# dmap$V1 = gsub("\\}", "", dmap$V1)
# 
# library(stringr)
# dmap = as.data.frame(str_split_fixed(dmap$V1, ",", 26))
# dmap = dmap[, c(2, 26, 1, 3, 9)]
# names(dmap) = c("chembl.id", "chembl.name", "phase", "indication", "ruleof5")
# dmap = dmap %>% mutate_all(as.character)
# dmap$chembl.name = gsub("[0-9]+.[0-9]+,", "", dmap$chembl.name)
# dmap[dmap == "None"] = NA
# save(dmap, file = file.path(dataFolder,"drug2715details.RData"))
# load(file.path(dataFolder,"drug2715details.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_output/spia_kegg_diseaseGenes.RData"))
load(file.path(dataFolder,"spia_output/spia_kegg_diseaseGenes_pseudo.RData"))

#'delete diseases which are not in both lists
Dis2Delete = setdiff(names(spia_kegg),names(spia_kegg_pseudo))
spia_kegg_tmp = spia_kegg

for (i in 1:length(spia_kegg)) {
  if(names(spia_kegg)[[i]] %in% Dis2Delete){
    spia_kegg_tmp[[i]] = NULL
  }
}

spia_kegg = spia_kegg_tmp
rm(Dis2Delete,spia_kegg_tmp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' prepare SPIA datadrame

spia_kegg = do.call(rbind, spia_kegg)
spia1 = data.table(spia_kegg[,c(5)])
names(spia1) = "pvalue"
spia1$spia = "Real"

spia_kegg_pseudo = do.call(rbind, spia_kegg_pseudo)
spia2= data.table(spia_kegg_pseudo[,c(5)])
names(spia2) = "pvalue"
spia2$spia = "Simulated"

spia = data.table(rbind(spia1,spia2))
rm(spia_kegg,spia_kegg_pseudo,spia1,spia2)
spia$spia = as.factor(spia$spia)
man_wtny <- wilcox.test(x = spia[spia == "Simulated", pvalue], y = spia[spia == "Real", pvalue])
print(man_wtny$p.value)
jpeg(file = file.path(dataFolder,"spia_real_pseudo.pvalues.boxplots.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
print(ggplot(spia, aes(x = spia, y = pvalue,fill=spia)) +
        geom_boxplot() +
        coord_cartesian(ylim = quantile(spia[, pvalue], c(0.03, 0.97))) +
        xlab("Distributions of the p-values from SPIA calculations") +
        ylab("p-value") +
        theme_bw(16) + theme(legend.position = "none")+
        scale_x_discrete(breaks = c("Simulated", "Real"), labels = c("Simulated Pathways", "KEGG Pathways")) +
        ggtitle(paste("p-value =", sprintf("%.2e", man_wtny$p.value))))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~ROC Curve~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# pred <- as.numeric(spia[, - log10(pvalue)])
# resp <- as.numeric(spia[, spia])
# 
# roc.curve <- roc(response = as.numeric(resp), predictor = as.numeric(pred), algorithm = 2, ci = TRUE, ci.method = "bootstrap", smooth = TRUE, boot.n = 100, parallel = TRUE, progress = "none")
# print(roc.curve)
# sp_ci <- ci.sp(roc.curve, sensitivities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
# se_ci <- ci.se(roc.curve, specifities = seq(0, 1, 0.05), boot.n = 1000, parallel = TRUE, progress = "none")
# jpeg(file = file.path(dataFolder,"disease_genes.roc.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
# par(pty = "s")
# plot(roc.curve, main = paste("AUC =", round(roc.curve$auc, 2)), xlab = "False positive rate", ylab = "True positive rate", identity.lty = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5)
# plot(se_ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
# plot(sp_ci, type = "shape", col = "lightgrey", border = NA, no.roc = TRUE)
# plot(roc.curve, add = TRUE, col = "#0066ff", lwd = 3)
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#