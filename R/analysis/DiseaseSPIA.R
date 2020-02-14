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
suppressWarnings(suppressMessages(library(VennDiagram)))
suppressWarnings(suppressMessages(library(RColorBrewer)))

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
#' Function to hide unnecessary outcome text in the terminal
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}
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

load(file.path(dataFolder,"geneID_v97.RData"))
gene_id$ENTREZ = gsub("^$", NA, gene_id$ENTREZ)
gene_id = gene_id[which(! is.na(gene_id$ENTREZ)),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lfc_comb = merge(lfc_comb, gene_id, by = "ensembl.id")
lfc_comb = lfc_comb[! duplicated(lfc_comb[, c('efo.id', 'ENTREZ')]),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' User can use any Disese response Dataset here with 5 columns:
#' "ensembl.id", "efo.id", "fold change", "ENTREZ", "HGNC"   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# lfc_efo = split(lfc_comb, lfc_comb$efo.id)
load(file.path(dataFolder,"results/disease_genes.RData"))
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

#' We like to discard diseases with less than 10 genes in the gene sets 
#' since after running SPIA or any other gene set enrichment method with 
#' such a small number of genes will not return almost any significant pathway.
#' However, users can comment out following line to not use the propsed filer
lfc_efo = Filter(function(x) dim(x)[1] > 10, lfc_efo) # remove diseases with very few (less than 10) genes to test

save(lfc_efo, file = file.path(dataFolder,"spia_input/lfc_disease_genes.RData"))
rm(lfc_comb)

#~~~~~~~~~~Create Named Vector for log fold changes in each disease~~~~~~~~~~#

# load(file.path(dataFolder,"spia_input/lfc_disease_genes.RData")) ###--get log fold change for all genes for each diseases-##

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

save(lfc_hgnc, lfc_ensembl, lfc_entrez, lfc_entrezID, hgnc_all, ensembl_all, entrez_all, entrezID_all, file = file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~Preapre SPIA Pathway Data sets~~~~~~~~~~~#

#~~~Prepare KEGG SPIA pathway file from KEGG pathway kgml files~~~#
#' Here we assume that the user obtained the KGML (xml) files from
#' KEGG ftp site (or downloaded them one by one from the KEGG website).
#' We will not keep KEGG pathway data since the access to the KEGG ftp server requires a license
#' However, KEGG license is open to use for academic purposes.

# makeSPIAdata(kgml.path=file.path(dataFolder,"spia_input/kgml/"),organism="hsa",out.path=file.path(dataFolder,"spia_input/real_kegg/"))

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~ KEGG SPIA ~~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))
rm(lfc_entrezID,lfc_hgnc,lfc_ensembl,ensembl_all,entrezID_all,hgnc_all)

pb <- txtProgressBar(min=0, max=length(lfc_entrez), style=3)
cat(sprintf("\n~~~~~SPIA calculation for Real KEGG Pathways~~~~~\n"))

spia_kegg = list()
for (i in 1 : length(lfc_entrez)) {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  spia_kegg[[i]] = quiet(spia(de = lfc_entrez[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa"))
}
close(pb)
names(spia_kegg) = names(lfc_entrez)

save(spia_kegg, file = file.path(dataFolder, "results/spia_output/spia_kegg_disease.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~ Reactome SPIA ~~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))
rm(lfc_entrez,lfc_hgnc,lfc_ensembl,ensembl_all,entrez_all,hgnc_all)

pb <- txtProgressBar(min=0, max=length(lfc_entrezID), style=3)
cat(sprintf("\n~~~~~SPIA calculation for Real Reactome Pathways~~~~~\n"))

spia_reactome = list()
for (i in 1 : length(lfc_entrezID)) {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  # cat(sprintf("\n~~~~~SPIA for Disease #%d~~~~~\n", (length(lfc_entrez) + 1) -i))
  spia_reactome[[i]] = quiet(spia(de = lfc_entrezID[[i]], all = entrezID_all, data.dir = file.path(dataFolder,"spia_input/real_react/"), organism = "hsa"))
}

close(pb)
names(spia_reactome) = names(lfc_entrezID)

save(spia_reactome, file = file.path(dataFolder, "results/spia_output/spia_reactome_disease.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~ Biocarta SPIA ~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))
rm(lfc_entrez,lfc_hgnc,lfc_ensembl,ensembl_all,entrez_all,hgnc_all)

pb <- txtProgressBar(min=0, max=length(lfc_entrezID), style=3)
cat(sprintf("\n~~~~~SPIA calculation for Real Biocarta Pathways~~~~~\n"))

spia_biocarta = list()
for (i in 1 : length(lfc_entrezID)) {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  # cat(sprintf("\n~~~~~SPIA for Disease #%d~~~~~\n", (length(lfc_entrez) + 1) -i))
  spia_biocarta[[i]] = quiet(spia(de = lfc_entrezID[[i]], all = entrezID_all, data.dir = file.path(dataFolder,"spia_input/real_biocarta/"), organism = "hsa"))
}

close(pb)
names(spia_biocarta) = names(lfc_entrezID)

save(spia_biocarta, file = file.path(dataFolder, "results/spia_output/spia_biocarta_disease.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~Simulated SPIA~~~~~~~~~~~#

#' we create simulated SPIA data sets for all three pathways
#' with random gene sets to see whether our results are by
#' random chance or meaningful indeed.


###___Create Simulated (Random) Pathways____####

#___create universe with ENTREZID:____#
load(file.path(dataFolder,"geneID_v97.RData"))
universe = unique(na.omit(gene_id$ENTREZ)) # for kegg
# universe = unique(na.omit(gsub("^","ENTREZID:",gene_id$ENTREZ))) # for reactome, biocarta

load(file.path(dataFolder,"spia_input/real_kegg/hsaSPIA.RData"))
# load(file.path(dataFolder,"spia_input/real_react/hsaSPIA.RData"))
# load(file.path(dataFolder,"spia_input/real_biocarta/hsaSPIA.RData"))

simulated_path = path.info
rm(path.info)

#_____Shuffle first 25 matrices in simulated_path______#
for (mylist in 1:length(simulated_path)) {
  for (mymat in 1:(length(simulated_path[[mylist]])-3)) { #not looping for last 3 indices since they are not matrices
    simulated_path[[mylist]][[mymat]][] = simulated_path[[mylist]][[mymat]][sample(length(simulated_path[[mylist]][[mymat]]))]
  }
}

#___Replace ENTREZ ID with ENTREZ ID from Universe_____#
for (mylist in 1:length(simulated_path)) {
  for (mymat in length(simulated_path[[mylist]])-2) { #need to change index number to '1' for reactome and biocarta, and to '2' for kegg
    simulated_path[[mylist]][[mymat]] = as.character(sample(universe,size = length(simulated_path[[mylist]][[mymat]]),replace = TRUE))
  }
}

#_____change row & column names of (25) matrices______#
for (mylist in 1:length(simulated_path)) {
  for (mymat in 1:(length(simulated_path[[mylist]])-3)) { #not looping for last 3 indices since they are not matrices
    colnames(simulated_path[[mylist]][[mymat]]) = simulated_path[[mylist]][[26]] #change index to '27' for reactome and biocarta, and to '26' for kegg
    rownames(simulated_path[[mylist]][[mymat]]) = simulated_path[[mylist]][[26]] #change index to '27' for reactome and biocarta, and to '26' for kegg
  }
}

path.info <- simulated_path

save(path.info,file = file.path(dataFolder,"spia_input/simulated_kegg/hsaSPIA.RData"))
# save(path.info,file = file.path(dataFolder,"spia_input/simulated_react/hsaSPIA.RData"))
# save(path.info,file = file.path(dataFolder,"spia_input/simulated_biocarta/hsaSPIA.RData"))
rm(path.info,simulated_path)

#_______________________________________________________##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~ Simulated KEGG SPIA ~~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))
rm(lfc_entrezID,lfc_hgnc,lfc_ensembl,ensembl_all,entrezID_all,hgnc_all)

pb <- txtProgressBar(min=0, max=length(lfc_entrez), style=3)
cat(sprintf("\n~~~~~SPIA calculation for simulated KEGG Pathways~~~~~\n"))

spia_kegg_simulated = list()
for (i in 1 : length(lfc_entrez)) {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  spia_kegg_simulated[[i]] = quiet(spia(de = lfc_entrez[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/simulated_kegg/"), organism = "hsa"))
}
close(pb)
names(spia_kegg_simulated) = names(lfc_entrez)

save(spia_kegg_simulated, file = file.path(dataFolder, "results/spia_output/spia_kegg_disease_simulated.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~ Simulated Reactome SPIA ~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))
rm(lfc_entrez,lfc_hgnc,lfc_ensembl,ensembl_all,entrez_all,hgnc_all)

pb <- txtProgressBar(min=0, max=length(lfc_entrezID), style=3)
cat(sprintf("\n~~~~~SPIA calculation for simulated Reactome Pathways~~~~~\n"))

spia_reactome_simulated = list()
for (i in 1 : length(lfc_entrezID)) {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  spia_reactome_simulated[[i]] = quiet(spia(de = lfc_entrezID[[i]], all = entrezID_all, data.dir = file.path(dataFolder,"spia_input/simulated_react/"), organism = "hsa"))
}
close(pb)
names(spia_reactome_simulated) = names(lfc_entrezID)

save(spia_reactome_simulated, file = file.path(dataFolder, "results/spia_output/spia_reactome_disease_simulated.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~ Simulated Biocarta SPIA ~~~~~~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder,"spia_input/lfc_disease_genes_namedVector.RData"))
rm(lfc_entrez,lfc_hgnc,lfc_ensembl,ensembl_all,entrez_all,hgnc_all)

pb <- txtProgressBar(min=0, max=length(lfc_entrezID), style=3)
cat(sprintf("\n~~~~~SPIA calculation for simulated Biocarta Pathways~~~~~\n"))

spia_biocarta_simulated = list()
for (i in 1 : length(lfc_entrezID)) {
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
  spia_biocarta_simulated[[i]] = quiet(spia(de = lfc_entrezID[[i]], all = entrezID_all, data.dir = file.path(dataFolder,"spia_input/simulated_biocarta/"), organism = "hsa"))
}
close(pb)
names(spia_biocarta_simulated) = names(lfc_entrezID)

save(spia_biocarta_simulated, file = file.path(dataFolder, "results/spia_output/spia_biocarta_disease_simulated.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~ Real_Vs_Simulated Pathway SPIA Boxplot Comparison ~~~~~#

#'~~~ Data for KEGG SPIA ~~~'#
load(file.path(dataFolder,"results/spia_output/spia_kegg_disease.RData"))
load(file.path(dataFolder,"results/spia_output/spia_kegg_disease_simulated.RData"))
spia = spia_kegg
spia_simulated = spia_kegg_simulated
rm(spia_kegg,spia_kegg_simulated)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#'~~~ Data for Reactome SPIA ~~~'#
# load(file.path(dataFolder,"results/spia_output/spia_reactome_disease.RData"))
# load(file.path(dataFolder,"results/spia_output/spia_reactome_disease_simulated.RData"))
# spia = spia_reactome
# spia_simulated = spia_reactome_simulated
# rm(spia_reactome,spia_reactome_simulated)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#'~~~ Data for Biocarta SPIA ~~~'#
# load(file.path(dataFolder,"results/spia_output/spia_biocarta_disease.RData"))
# load(file.path(dataFolder,"results/spia_output/spia_biocarta_disease_simulated.RData"))
# spia = spia_biocarta
# spia_simulated = spia_biocarta_simulated
# rm(spia_biocarta,spia_biocarta_simulated)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#'delete diseases which are not in both lists
Dis2Delete = setdiff(names(spia),names(spia_simulated))
spia_tmp = spia

for (i in 1:length(spia)) {
  if(names(spia)[[i]] %in% Dis2Delete){
    spia_tmp[[i]] = NULL
  }
}

spia = spia_tmp
rm(Dis2Delete,spia_tmp)

#' prepare SPIA dataframe
spia = do.call(rbind, spia)
spia1 = data.table(spia[,c(5)])
names(spia1) = "pvalue"
spia1$spia = "Real"

spia_simulated = do.call(rbind, spia_simulated)
spia2= data.table(spia_simulated[,c(5)])
names(spia2) = "pvalue"
spia2$spia = "Simulated"

spia_full = data.table(rbind(spia1,spia2))
rm(spia,spia_simulated,spia1,spia2)
spia_full$spia = as.factor(spia_full$spia)
man_wtny <- wilcox.test(x = spia_full[spia == "Simulated", pvalue], y = spia_full[spia == "Real", pvalue])
print(man_wtny$p.value)
jpeg(file = file.path(dataFolder, "results/figures/spia_real_simulated_pvalues_violinplots_kegg.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
# jpeg(file = file.path(dataFolder, "results/figures/spia_real_simulated_pvalues_violinplots_reactome.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
# jpeg(file = file.path(dataFolder, "results/figures/spia_real_simulated_pvalues_violinplots_biocarta.jpeg"), width = 6 * 150, height = 6 * 150, res = 150)
ggplot(spia_full, aes(x = spia, y = pvalue,fill=spia)) +
        geom_violin(trim=FALSE) + geom_boxplot(width=0.1)+
        xlab("Distributions of the p-values from SPIA calculations") +
        ylab("p-value") +
        theme_bw(16) + theme(legend.position = "none")+
        scale_x_discrete(breaks = c("Simulated", "Real"), labels = c("Simulated KEGG Pathways", "KEGG Pathways")) +
        # scale_x_discrete(breaks = c("Simulated", "Real"), labels = c("Simulated Reactome Pathways", "Reactome Pathways")) +
        # scale_x_discrete(breaks = c("Simulated", "Real"), labels = c("Simulated Biocarta Pathways", "Biocarta Pathways")) +
        ggtitle(paste("p-value =", sprintf("%.2e", man_wtny$p.value)))

dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~: Pairwise Gene sets intersection :~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Compute all pairwise intersections between any given 'Disease Gene sets'.
pairwise_genesets <- list()

for (i in 1:(length(lfc_efo)-1)) {
  ith_intersection <- list()
  ith_percentile <- list()
  for (j in (i+1):length(lfc_efo)) {
    ith_intersection[[ names(lfc_efo)[j] ]] <- intersect(lfc_efo[[i]]$ensembl.id, lfc_efo[[j]]$ensembl.id)
    # ith_percentile[[ names(lfc_efo)[j] ]] <- round((length(intersect(lfc_efo[[i]]$ensembl.id, lfc_efo[[j]]$ensembl.id))/length(lfc_efo[[i]]$ensembl.id))*100,2)
  }
  pairwise_genesets[[ names(lfc_efo)[i] ]] <- ith_intersection
  # pairwise_genesets[[ names(lfc_efo)[i] ]] <- ith_percentile
}

save(pairwise_genesets, file = file.path(dataFolder, "results/spia_output/pairwise_genesets_intersection.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~: Disease Gene sets intersection with Pathway Gene sets:~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~: Disease Gene Sets :~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Disease Gene lists
load(file.path(dataFolder,"spia_input/lfc_disease_genes.RData"))

DisGeneLists = list()
for (i in seq_along(lfc_efo)) {
  DisGeneLists[[i]] = unlist(lfc_efo[[i]][["ENTREZ"]],recursive = F)
  #names(DisGeneLists)[[i]] = names(lfc_efo)[[i]]
}

# Disease Gene combined list
DisGeneList = unique(unlist(DisGeneLists, recursive = F))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~: KEGG Gene Sets :~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Pathway Gene lists
load(file.path(dataFolder,"spia_input/real_kegg/hsaSPIA.RData"))
PwGeneLists = list()
for (i in seq_along(path.info)) {
  PwGeneLists[[i]] = unlist(path.info[[i]][["nodes"]],recursive = F)
  #names(PwGeneLists)[[i]] = path.info[[i]][["title"]]
}
rm(path.info)
# Pathway Gene combined list
PwGeneListKegg = unique(unlist(PwGeneLists, recursive = F)) # 5536 in KEGG, 10000 in Reactome,

# Calculate Disease Genes that are not in any of KEGG patways
DisGenesNotInanyPW = setdiff(DisGeneList,PwGeneListKegg) # 2031 in KEGG
x = intersect(DisGeneList,PwGeneListKegg)
cat(sprintf("'Number of Disease Genes which are not present in any of the KEGG pathways : %d",length(DisGenesNotInanyPW)))
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~: Reactome Gene Sets :~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Pathway Gene lists
load(file.path(dataFolder,"spia_input/real_react/hsaSPIA.RData"))
PwGeneLists = list()
for (i in seq_along(path.info)) {
  PwGeneLists[[i]] = unlist(path.info[[i]][["nodes"]],recursive = F)
  #names(PwGeneLists)[[i]] = path.info[[i]][["title"]]
}
rm(path.info)
# Pathway Gene combined list
PwGeneListReact = unique(unlist(PwGeneLists, recursive = F)) # 5536 in KEGG, 10000 in Reactome,
PwGeneListReact = gsub("^ENTREZID:","",PwGeneListReact) # for cleaning reactome and Biocarta IDs

# Calculate Disease Genes that are not in any of KEGG patways
DisGenesNotInanyPW = setdiff(DisGeneList,PwGeneListReact) # 1267 in Reactome
x = intersect(DisGeneList,PwGeneListReact)
cat(sprintf("'Number of Disease Genes which are not present in any of the Reactome pathways : %d",length(DisGenesNotInanyPW)))
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~: Biocarta Gene Sets :~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Pathway Gene lists
load(file.path(dataFolder,"spia_input/real_biocarta/hsaSPIA.RData"))
PwGeneLists = list()
for (i in seq_along(path.info)) {
  PwGeneLists[[i]] = unlist(path.info[[i]][["nodes"]],recursive = F)
  #names(PwGeneLists)[[i]] = path.info[[i]][["title"]]
}
rm(path.info)
# Pathway Gene combined list
PwGeneListBiocarta = unique(unlist(PwGeneLists, recursive = F)) # 5536 in KEGG, 10000 in Reactome,
PwGeneListBiocarta = gsub("^ENTREZID:","",PwGeneListBiocarta) # for cleaning reactome and Biocarta IDs

# Calculate Disease Genes that are not in any of KEGG patways
DisGenesNotInanyPW = setdiff(DisGeneList,PwGeneListBiocarta) # 2895 in BIocarta
x = intersect(DisGeneList,PwGeneListBiocarta)
cat(sprintf("'Number of Disease Genes which are not present in any of the Biocarta pathways : %d",length(DisGenesNotInanyPW)))
gc()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate Disease Genes that are present in more than five of KEGG patways
DisGenFreq = as.list(rep(0, length(DisGeneList)))

for (i in seq_along(DisGeneList)) {
  for (j in seq_along(PwGeneLists)) {
    PwGeneLists[[j]] = gsub("^ENTREZID:","",PwGeneLists[[j]]) # clean gene IDs before comparison in Reactome and Biocarta
    if (DisGeneList[[i]] %in% PwGeneLists[[j]]) {
      DisGenFreq[[i]] = DisGenFreq[[i]]+1
    }
  }
  names(DisGenFreq)[[i]] = as.character(DisGeneList[[i]])
}
sum(DisGenFreq > 5) # 272 in KEGG, 940 in Reactome and 36 in Biocarta
cat(sprintf("'Number of Disease Genes which are present in more than 5 pathways : %d",sum(DisGenFreq > 5)))
gc()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

save(DisGeneList,PwGeneListKegg,PwGeneListReact,PwGeneListBiocarta, file=file.path(dataFolder, "results/DisPwGenList.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Venn Diagram for Gene Intersection among the Disease and Pathway Genesets
load(file.path(dataFolder, "results/DisPwGeneList.RData"))

myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(DisGeneList,PwGeneListKegg, PwGeneListReact,PwGeneListBiocarta),
  category.names = c("Disease Geneset","KEGG Geneset", "Reactome Geneset","Biocarta Geneset"),
  filename = file.path(dataFolder, 'results/figures/Geneset_intersection_venn_diagramm.png'),
  output=TRUE,
  main = "Gene Intersections among the Disease and Pathway Genesets",
  main.cex = 2,
  
  # Output features
  imagetype="png" ,
  height = 3000 , 
  width = 3000 , 
  resolution = 250,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
