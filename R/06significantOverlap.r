#' 6th script
#' summary:
#' calculate significant overlap between/among different data sets 
#' e.g., Drug Perturbed Genes & DEGs data,
#' Drug Perturbed Genes & GWAS data,
#' GWAS & DEGs data,
#' Drug Perturbed Genes & DEGs & GWAS data

library(EnsDb.Hsapiens.v86)

library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

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
opentargets.drugs <- unique(fread("./data/opentargets.drugs.tsv"))
opentargets.degs <- unique(fread("./data/opentargets.degs.tsv"))

#################################################


#-------------------------------------------------#
##____DEGs to GWAS Genes >>> DISEASE Genes_____####
#_________________________________________________#


# split by disease
stopgap.list <- split(stopgap, stopgap$efo.id)
opentargets.degs.list <- split(opentargets.degs, opentargets.degs$efo.id)

# create universes for Fisher's test
ensembl.ids <- unique(keys(EnsDb.Hsapiens.v86))
efo.ids <- unique(c(stopgap[, efo.id], opentargets.degs[, efo.id], opentargets.drugs[, efo.id]))
chembl.ids <- unique(opentargets.drugs[, chembl.id])

#### Get Gene Symbol for Ensembl IDs 
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", version=93, dataset="hsapiens_gene_ensembl")
val <- c(1:22,'X','Y','MT')
chr_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name'), 
                   filters ='chromosome_name', values =val, mart = ensembl)
chr_genes$hgnc_symbol <- gsub("^$", NA, chr_genes$hgnc_symbol)
chr_genes <- chr_genes[which(!is.na(chr_genes$hgnc_symbol)),]
chr_genes <- unique(chr_genes)
hgnc.list = chr_genes[[1]]

setdiff(ensembl.ids,hgnc.list)

## Fisher's test
# loop through diseases

fisher.genes <- foreach (i = seq(opentargets.degs.list), .combine = rbind, .errorhandling = "remove") %do% {
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

fisher.genes <- unique(fisher.genes)

# correct p-values
fisher.genes[p.value == 0, p.value := 3e-324]
fisher.genes <- fisher.genes[order(p.value), ]
fisher.genes[, p.adjusted := p.adjust(p.value, method = "fdr")]

# add dummy variable
fisher.genes[, same.disease := ifelse(efo.id.opentargets.degs == efo.id.stopgap, TRUE, FALSE)]
fwrite(fisher.genes, "./data/fisher.genes.commongenes.tsv", sep = "\t")


##----------------------------------##
##____Drugs to DISEASE Genes_____#####
##----------------------------------##
opentargets.drugs <- unique(fread("./data/opentargets.drugs.tsv"))
harmonizome <- unique(fread("./data/harmonizome.tsv"))
chembl.ids <- unique(harmonizome[, chembl.id])

##### create Disease-Genes list #####
fisher.genes <- fread("./data/fisher.genes.commongenes.tsv")
DisGen = fisher.genes[same.disease == TRUE & overlap > 0]
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.opentargets.degs),] #remove duplicated rows based on one column
DisGenX <- unique(DisGen %>% 
                         mutate(commonGenes = strsplit(as.character(commonGenes), "\\|")) %>% 
                         unnest(commonGenes))
DisGen.list <- split(DisGenX, DisGenX$efo.id.opentargets.degs)
efo.ids <- unique(DisGen$efo.id.opentargets.degs)

###___Fisher's test____#####
super.drugs <- foreach (i = seq(efo.ids), .combine = rbind, .errorhandling = "remove") %dopar% {
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
      #tmp <- merge(unique(opentargets.drugs[, .(chembl.id, chembl.name)]), tmp, by = "chembl.id")
      # add dummy variable
      tmp[, existing.indication := ifelse(nrow(opentargets.drugs[efo.id == this.efo.id & chembl.id == this.chembl.id]) > 0, TRUE, FALSE)]
      setcolorder(tmp, c("chembl.id", "efo.id", "efo.term", "DEGs", "GWAGs", "overlap", "universe","commonGenes", "odds.ratio", "p.value", "existing.indication"))
    }
  }
}


super.drugs <- unique(super.drugs)
super.drugs = super.drugs[overlap > 0]
length(unique(super.drugs$chembl.id))
length(unique(super.drugs$commonGenes))
#save(super.drugs,file="./data/super.drugs.48D.RData")

load(("./data/super.drugs.48D.RData"))
# correct p-values
super.drugs[, p.adjusted := p.adjust(p.value, method = "fdr")]
super.drugs <- super.drugs[p.adjusted < 0.05]
#save(super.drugs,file="./data/super.drugs.48D.padj.RData")
super.drugs <- super.drugs[p.adjusted < 1e-05]
save(super.drugs,file="./data/super.drugs.48D.padj1e-5.RData")
#md2 <- unique(min.drugs[,c(1,3,12)]) #filtering columns for merging indication area to our super drugs
#super.drugs <- merge(super.drugs,md2,by=c("chembl.id","efo.id"))
super.drugs <- super.drugs[order(p.adjusted), ]





# #__________________________#
#####___Drug Data Preperation___#####

load("./data/gene.id.entrez.RData")

# fisher.drugs = fread("./data/fisher.drugs.commongenes.tsv")
# fisher.drugs = as.data.frame(unique(fisher.drugs$chembl.id))
# names(fisher.drugs) = "chembl.id"

harmonizome = unique(fread("./data/harmonizome.tsv"))
harmonizome = harmonizome[,c(5,1,7)]
length(unique(harmonizome$chembl.id))
length(unique(harmonizome$ensembl.id))
harmonizome = harmonizome[order(ensembl.id,decreasing = TRUE), ]
harmonizome = harmonizome[!duplicated(harmonizome[,c('ensembl.id', 'chembl.id')]),]
harmonizome = merge(harmonizome,gene.id.entrez,by="ensembl.id")
#harmonizome = merge(fisher.drugs,harmonizome,by="chembl.id")
# harmonizome = harmonizome[,c(1,2,4,3)]
#harmonizome.chembl = split(harmonizome, harmonizome$chembl.id)

# load("./data/super.drugs.RData")
# super.drugs = super.drugs[,c(1,2,9)]
load("./data/super.drugs.48D.padj.RData")
load("./data/super.drugs.48D.padj1e-5.RData")
super.drugs = super.drugs[,c(1,2,8)]

#super.drugs.chembl = split(super.drugs, super.drugs$chembl.id)

drug.genes <- unique(super.drugs %>% 
                       mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>% 
                       unnest(ensembl.id))
drug.genes$ensembl.id = gsub("\"","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("c\\(","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("\\)","",drug.genes$ensembl.id)
drug.genes$ensembl.id = trimws(drug.genes$ensembl.id)

length(unique(drug.genes$chembl.id))
length(unique(drug.genes$ensembl.id))
drug.genes = merge(drug.genes,harmonizome,by=c('chembl.id','ensembl.id')) # merge with harmonizome
drug.genes = drug.genes[!duplicated(drug.genes[,c('chembl.id','ensembl.id','efo.id')]),]

drug.genes.list <- split(drug.genes, list(drug.genes$chembl.id,drug.genes$efo.id)) #
drug.genes.list= Filter(function(x) dim(x)[1] > 0, drug.genes.list) # remove empty lists
save(drug.genes.list,file = "./data/drug.genes.list.48D.padj.RData")
save(drug.genes.list,file = "./data/drug.genes.list.48D.padj1e-5.RData")

####____Create named entity list for each drug_disease pairs___####
# load("./data/drug.genes.list.RData")
load("./data/drug.genes.list.48D.padj.RData")
load("./data/drug.genes.list.48D.padj1e-5.RData")

lfc_drug_ensembl <- list()
for (i in 1:length(drug.genes.list)) {
  lfc_drug_ensembl[[i]] = setNames(as.numeric(drug.genes.list[[i]][[4]]),as.character(drug.genes.list[[i]][[2]]))
  names(lfc_drug_ensembl)[i] = names(drug.genes.list)[i]
}

lfc_drug_entrez <- list()
for (i in 1:length(drug.genes.list)) {
  lfc_drug_entrez[[i]] = setNames(as.numeric(drug.genes.list[[i]][[5]]),as.character(drug.genes.list[[i]][[4]]))
  names(lfc_drug_entrez)[i] = names(drug.genes.list)[i]
}

lfc_drug_entrezID <- list()
for (i in 1:length(drug.genes.list)) {
  lfc_drug_entrezID[[i]] = setNames(as.numeric(drug.genes.list[[i]][[5]]),as.character(gsub("^","ENTREZID:",drug.genes.list[[i]][[4]])))
  names(lfc_drug_entrezID)[i] = names(drug.genes.list)[i]
}

##_____Create a vector with all Gene universe to proxy Array Genes_________###
load("./data/gene.id.entrez.RData")
ensembl_all = unique(gene.id.entrez$ensembl.id)
entrez_all = unique(gene.id.entrez$ENTREZ)
entrezID_all = unique(gsub("^","ENTREZID:",gene.id.entrez$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_drug_ensembl,ensembl_all,file="./data/lfc.drug.48D.padj.ensembl.RData")
save(lfc_drug_entrez,entrez_all,file="./data/lfc.drug.48D.padj.entrez.RData")
save(lfc_drug_entrezID,entrezID_all,file="./data/lfc.drug.48D.padj.entrezID.RData")

save(lfc_drug_ensembl,ensembl_all,file="./data/lfc.drug.48D.padj1e-5.ensembl.RData")
save(lfc_drug_entrez,entrez_all,file="./data/lfc.drug.48D.padj1e-5.entrez.RData")
save(lfc_drug_entrezID,entrezID_all,file="./data/lfc.drug.48D.padj1e-5.entrezID.RData")
##____________________________________________________##




##-----------------------------##
#####_____Pathway SPIA______#####
##-----------------------------##


#####__Prepare DEG data set for SPIA__#####
degs <- unique(fread("./data/opentargets.degs.tsv"))
degs <- degs[,c(1,3,6,7)]

# process positive lfc 
lfc.pos <- degs[lfc >= 0]
#lfc.pos <- lfc.pos[order(lfc,decreasing = TRUE), ] #lfc or pvalue we should sort the dataframe?
lfc.pos$pval.pos <- abs(log10(lfc.pos$pval)) #created dummy variable for pvalue to order it decreasingly, since there are multiple lfc with same pvalue sometimes.
#lfc.pos <- lfc.pos[order(pval), ]
lfc.pos <- lfc.pos[order(lfc,pval.pos,decreasing = TRUE), ] 
lfc.pos = lfc.pos[!duplicated(lfc.pos[,c('efo.id', 'ensembl.id')]),]
#lfc.pos.efo <- split(lfc.pos, lfc.pos$efo.id)
lfc.pos$pval.pos <- NULL

# process negative lfc 
lfc.neg = degs[lfc < 0]
#lfc.neg = lfc.neg[order(lfc), ] #lfc or pvalue we should sort the dataframe?
lfc.neg = lfc.neg[order(pval,lfc), ]
lfc.neg = lfc.neg[!duplicated(lfc.neg[,c('efo.id','ensembl.id')]),]
#lfc.neg.efo = split(lfc.neg, lfc.neg$efo.id)

# combine both lfc 
lfc.com = rbind(lfc.pos,lfc.neg) #combine positive and negative lfc change back to a single data frame
lfc.com = lfc.com[order(pval), ]
lfc.com = lfc.com[!duplicated(lfc.com[,c('efo.id','ensembl.id')]),]
lfc.com$pval <- NULL
rm(lfc.neg,lfc.pos,degs)
#__map ensembl IDs to ENTREZ ID___#

#-----get mapping among ENTREZ_HGNC_ENSEMBL_IDs---------#
# library(biomaRt)
# ensembl = useEnsembl(biomart="ensembl", version=92, dataset="hsapiens_gene_ensembl") # v92 has less id than v79
# val <- c(1:23,"X","Y")
# gene.id <- getBM(attributes=c('entrezgene','hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position'),
#                    filters ='chromosome_name', values =val, mart = ensembl)
# 
# save(gene.id,file = "./data/geneID.RData")
# rm(ensembl,val)

load("./data/geneID.RData")
gene.id$entrezgene <- gsub("^$", NA, gene.id$entrezgene)
gene.id <- gene.id[which(!is.na(gene.id$entrezgene)),]
gene.id <- data.table(gene.id)
gene.id<- unique(gene.id[,c('entrezgene','ensembl_gene_id')])
names(gene.id) <- c("ENTREZ","ensembl.id")
gene.id <- gene.id[!duplicated(gene.id$ensembl.id),]

lfc.com <- merge(lfc.com,gene.id,by="ensembl.id")
lfc.com = lfc.com[!duplicated(lfc.com[,c('efo.id','ENTREZ')]),]



fisher.genes <- fread("./data/fisher.genes.commongenes.tsv")
DisGen = fisher.genes[same.disease == TRUE & overlap > 0]
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.opentargets.degs),] #remove duplicated rows based on one column
DisGen = DisGen[,c(1,2,9)]
DisGen <- unique(DisGen %>% 
                       mutate(ensembl.id = strsplit(as.character(commonGenes), "\\|")) %>% 
                       unnest(ensembl.id))
DisGen$ensembl.id = gsub("\"","",DisGen$ensembl.id)
DisGen$ensembl.id = gsub("c\\(","",DisGen$ensembl.id)
DisGen$ensembl.id = gsub("\\)","",DisGen$ensembl.id)
DisGen$ensembl.id = trimws(DisGen$ensembl.id)
length(unique(DisGen$chembl.id))
length(unique(DisGen$ensembl.id))
DisGen$commonGenes = NULL
names(DisGen) = c("efo.id","efo.term","ensembl.id")
DisGen = merge(DisGen,lfc.com,by=c('ensembl.id','efo.id')) # merge with harmonizome

lfc.efo = split(DisGen, DisGen$efo.id)
rm(degs,lfc.neg,lfc.pos,lfc.com)

save(lfc.efo,file="./data/lfc.dis52.RData")

##_____read log fold changes for a particular disease_______###
load("./data/lfc.dis52.RData") ###--get log fold change for all genes for each diseases-##

##_____Create Named Vector for log fold changes in each disease_____________###

lfc_ensembl <- list()
for (i in 1:length(lfc.efo)) {
  lfc_ensembl[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]),as.character(lfc.efo[[i]][[1]]))
  names(lfc_ensembl)[[i]] = names(lfc.efo)[[i]]
}

lfc_entrez <- list()
for (i in 1:length(lfc.efo)) {
  lfc_entrez[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]),as.character(lfc.efo[[i]][[5]]))
  names(lfc_entrez)[[i]] = names(lfc.efo)[[i]]
}

lfc_entrezID <- list()
for (i in 1:length(lfc.efo)) {
  lfc_entrezID[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]),as.character(gsub("^","ENTREZID:",lfc.efo[[i]][[5]])))
  names(lfc_entrezID)[[i]] = names(lfc.efo)[[i]]
}


##_____Create a vector with all Gene universe to proxy Array Genes_________###
ensembl_all = unique(gene.id$ensembl.id)
entrez_all = unique(gene.id$ENTREZ)
entrezID_all = unique(gsub("^","ENTREZID:",gene.id$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_ensembl,lfc_entrez,lfc_entrezID,ensembl_all,entrez_all,entrezID_all,file="./data/lfc.namedVec.dis52.RData")
rm(lfc.efo,gene.id)

#################################3

#####________SPIA____________#####

library(SPIA)
#___________KEGG SPIA______________#

load("./data/lfc.namedVec.dis52.RData")
spia_kegg = list()
for (i in 1:length(lfc_entrez)) {
  spia_kegg[[i]] = spia(de = lfc_entrez[[i]], all = entrez_all, data.dir="./data/real_kegg/",organism="hsa")
}
names(spia_kegg) = names(lfc.efo)
save(spia_kegg,file = "./data/spia/spia_kegg_dis52_results.RData")
# rm(list=ls())
# gc()

#load("./data/lfc.all.RData")
spia_kegg_fake = list()
for (i in 1:length(lfc_entrez)) {
  spia_kegg_fake[[i]] = spia(de = lfc_entrez[[i]], all = entrez_all, data.dir="./data/pseudo_kegg/",organism="hsa")
}
save(spia_kegg_fake,file = "./data/spia/spia_kegg_dis52_fake_results.RData")
plotP(spia_kegg[[13]])
rm(list=ls())
gc()

load("./data/spia/spia_kegg_dis52_results.RData")

spia_kegg_ad = spia(de = lfc_entrez[[3]], all = entrez_all, data.dir="./data/real_kegg/",organism="hsa")
plotP(spia_kegg_ad)

# pdf(file="./data/spia/spia_kegg_plotsX.pdf")
# plotP(spia_kegg,threshold=0.1)
# plotP(spia_kegg_fake,threshold=0.1)
# dev.off()
# 
# spia_kegg1 <- data.table(spia_kegg)
# spia_kegg1 <- spia_kegg1[pGFdr <= 0.05]
# spia_kegg_fake1 <- data.table(spia_kegg_fake)
# spia_kegg_fake1 <- spia_kegg_fake1[pGFdr <= 0.05] 
##################################################################



