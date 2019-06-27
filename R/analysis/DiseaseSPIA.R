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
lfc.pos = DEGs[lfc >= 0]
#lfc.pos = lfc.pos[order(lfc,decreasing = TRUE), ] #lfc or pvalue we should sort the dataframe?
lfc.pos$pval.pos = abs(log10(lfc.pos$pval)) #created dummy variable for pvalue to order it decreasingly, since there are multiple lfc with same pvalue sometimes.
#lfc.pos = lfc.pos[order(pval), ]
lfc.pos = lfc.pos[order(lfc, pval.pos, decreasing = TRUE),]
lfc.pos = lfc.pos[! duplicated(lfc.pos[, c('efo.id', 'ensembl.id')]),]
#lfc.pos.efo = split(lfc.pos, lfc.pos$efo.id)
lfc.pos$pval.pos = NULL

# process negative lfc
lfc.neg = DEGs[lfc < 0]
#lfc.neg = lfc.neg[order(lfc), ] #lfc or pvalue we should sort the dataframe?
lfc.neg = lfc.neg[order(pval, lfc),]
lfc.neg = lfc.neg[! duplicated(lfc.neg[, c('efo.id', 'ensembl.id')]),]
#lfc.neg.efo = split(lfc.neg, lfc.neg$efo.id)

# combine both lfc
lfc.com = rbind(lfc.pos, lfc.neg) #combine positive and negative lfc change back to a single data frame
lfc.com = lfc.com[order(pval),]
lfc.com = lfc.com[! duplicated(lfc.com[, c('efo.id', 'ensembl.id')]),]
lfc.com$pval = NULL
rm(lfc.neg, lfc.pos, DEGs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"geneID.RData"))
gene.id$entrezgene = gsub("^$", NA, gene.id$entrezgene)
gene.id = gene.id[which(! is.na(gene.id$entrezgene)),]
gene.id = data.table(gene.id)
gene.id = unique(gene.id[, c('entrezgene', 'ensembl_gene_id', 'hgnc_symbol')])
names(gene.id) = c("ENTREZ", "ensembl.id", "HGNC")
gene.id = gene.id[! duplicated(gene.id$ensembl.id),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lfc.com = merge(lfc.com, gene.id, by = "ensembl.id")
lfc.com = lfc.com[! duplicated(lfc.com[, c('efo.id', 'ENTREZ')]),]

# lfc.efo = split(lfc.com, lfc.com$efo.id)
load(file.path(dataFolder,"disease.genes50.RData"))
DisGen = disease.genes[same.disease == TRUE &
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
length(unique(DisGen$efo.id))
length(unique(DisGen$ensembl.id))
DisGen$commonGenes = NULL
names(DisGen) = c("efo.id", "efo.term", "ensembl.id")
DisGen = merge(DisGen, lfc.com, by = c('ensembl.id', 'efo.id')) # merge with harmonizome

lfc.efo = split(DisGen, DisGen$efo.term)
lfc.efo = Filter(function(x) dim(x)[1] > 10, lfc.efo) # remove diseases with very few (less than 10) genes to test

save(lfc.efo, file = file.path(dataFolder,"disease47.genes50.lfc.RData"))
rm(lfc.com)

#~~~~~~~~~~Create Named Vector for log fold changes in each disease~~~~~~~~~~#

load(file.path(dataFolder,"disease47.genes50.lfc.RData"))###--get log fold change for all genes for each diseases-##

lfc_ensembl = list()
for (i in 1 : length(lfc.efo)) {
    lfc_ensembl[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]), as.character(lfc.efo[[i]][[1]]))
    names(lfc_ensembl)[[i]] = names(lfc.efo)[[i]]
}

lfc_entrez = list()
for (i in 1 : length(lfc.efo)) {
    lfc_entrez[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]), as.character(lfc.efo[[i]][[5]]))
    names(lfc_entrez)[[i]] = names(lfc.efo)[[i]]
}

lfc_entrezID = list()
for (i in 1 : length(lfc.efo)) {
    lfc_entrezID[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]), as.character(gsub("^", "ENTREZID:", lfc.efo[[i]][[5]])))
    names(lfc_entrezID)[[i]] = names(lfc.efo)[[i]]
}

lfc_hgnc = list()
for (i in 1 : length(lfc.efo)) {
    lfc_hgnc[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]), as.character(lfc.efo[[i]][[6]]))
    names(lfc_hgnc)[[i]] = names(lfc.efo)[[i]]
}

#~~~~~~Create a vector with all Gene universe to proxy Array Genes~~~~~~~#
hgnc_all = unique(gene.id$HGNC)
ensembl_all = unique(gene.id$ensembl.id)
entrez_all = unique(gene.id$ENTREZ)
entrezID_all = unique(gsub("^", "ENTREZID:", gene.id$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_hgnc, lfc_ensembl, lfc_entrez, lfc_entrezID, hgnc_all, ensembl_all, entrez_all, entrezID_all, file = file.path(dataFolder,"disease47.genes50.lfc.namedVec.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~Preapre SPIA Pathway Data sets~~~~~~~~~~~#

#~~~Prepare KEGG SPIA pathway file from KEGG pathway kgml files~~~#
#' Here we assume that the user obtained the KGML (xml) files from
#' KEGG ftp site (or downloaded them one by one from the KEGG website).
#' We will not keep KEGG pathway data since the access to the KEGG ftp server requires a license

makeSPIAdata(kgml.path=file.path(dataFolder,"spia_input/kgml/"),organism="hsa",out.path=file.path(dataFolder,"spia_input/real_kegg/"))

#~~~Prepare Reactome SPIA pathway file with graphite package~~~#

reactome <- pathways("hsapiens", "reactome")
reactome <- convertIdentifiers(reactome, "ENTREZID")
prepareSPIA(reactome, file.path(dataFolder,"spia_input/real_react/hsa"))

#~~~Prepare Reactome SPIA pathway file with graphite package~~~#
biocarta <- pathways("hsapiens", "biocarta")
biocarta <- convertIdentifiers(biocarta, "ENTREZID")
prepareSPIA(biocarta, file.path(dataFolder,"spia_input/real_biocarta/hsa"))

##_______________________________###


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~ Disease SPIA~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' SPIA can be calcualted for KEGG, Reactome and WikiPathways
#' Here, we showed only with KEGG Pathways
#' Necessary data sets for calculating SPIA for Reactome
#' and WikiPathways are uploaded in the git repository


#~~~~~~~~~~~~~KEGG SPIA~~~~~~~~~~~~~#

load(file.path(dataFolder,"disease47.genes50.lfc.namedVec.RData"))

spia_kegg = list()
for (i in 1 : length(lfc_entrez)) {
    spia_kegg[[i]] = spia(de = lfc_entrez[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}
names(spia_kegg) = names(lfc_entrez)

save(spia_kegg, file = file.path(dataFolder,"spia_output/spia_kegg_disease47.genes50_results.RData"))

spia_kegg_degs = lapply(spia_kegg_degs, function(x) x[x$pNDE <= 0.05,])
plotP(spia_kegg_degs[[4]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~Pseudo KEGG SPIA~~~~~~~~~~~#

#' we create pseudo SPIA data sets for all three pathways
#' with random gene sets to see whether our results are by
#' random chance or meaningful indeed.


###___Create Pseudo (Random) Pathways____####

#___create universe with ENTREZID:____#
load(file.path(dataFolder,"geneID.RData"))
universe = unique(na.omit(gene.id$entrezgene)) # for kegg
# universe = unique(na.omit(gsub("^","ENTREZID:",gene.id$entrezgene))) # for reactome, biocarta

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

spia_kegg_pseudo = list()
for (i in 1 : length(lfc_entrez)) {
    spia_kegg_pseudo[[i]] = spia(de = lfc_entrez[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/pseudo_kegg/"), organism = "hsa")
}
save(spia_kegg_pseudo, file = file.path(dataFolder,"spia_output/spia_kegg_disease47.genes50_pseudo_results.RData"))
rm(list = ls())
gc()


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


