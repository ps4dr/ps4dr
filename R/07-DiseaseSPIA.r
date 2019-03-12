library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

setwd("/home/memon/projects/msdrp/")

##-----------------------------##
#####_____Pathway SPIA______#####
##-----------------------------##


#####__Prepare DEG data set for SPIA__#####
load("./data/DEGs.RData")
DEGs <- DEGs[,c(1,3,6,7)]

# process positive lfc 
lfc.pos <- DEGs[lfc >= 0]
#lfc.pos <- lfc.pos[order(lfc,decreasing = TRUE), ] #lfc or pvalue we should sort the dataframe?
lfc.pos$pval.pos <- abs(log10(lfc.pos$pval)) #created dummy variable for pvalue to order it decreasingly, since there are multiple lfc with same pvalue sometimes.
#lfc.pos <- lfc.pos[order(pval), ]
lfc.pos <- lfc.pos[order(lfc,pval.pos,decreasing = TRUE), ] 
lfc.pos = lfc.pos[!duplicated(lfc.pos[,c('efo.id', 'ensembl.id')]),]
#lfc.pos.efo <- split(lfc.pos, lfc.pos$efo.id)
lfc.pos$pval.pos <- NULL

# process negative lfc 
lfc.neg = DEGs[lfc < 0]
#lfc.neg = lfc.neg[order(lfc), ] #lfc or pvalue we should sort the dataframe?
lfc.neg = lfc.neg[order(pval,lfc), ]
lfc.neg = lfc.neg[!duplicated(lfc.neg[,c('efo.id','ensembl.id')]),]
#lfc.neg.efo = split(lfc.neg, lfc.neg$efo.id)

# combine both lfc 
lfc.com = rbind(lfc.pos,lfc.neg) #combine positive and negative lfc change back to a single data frame
lfc.com = lfc.com[order(pval), ]
lfc.com = lfc.com[!duplicated(lfc.com[,c('efo.id','ensembl.id')]),]
lfc.com$pval <- NULL
rm(lfc.neg,lfc.pos,DEGs)
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
gene.id<- unique(gene.id[,c('entrezgene','ensembl_gene_id','hgnc_symbol')])
names(gene.id) <- c("ENTREZ","ensembl.id","HGNC")
gene.id <- gene.id[!duplicated(gene.id$ensembl.id),]


lfc.com <- merge(lfc.com,gene.id,by="ensembl.id")
lfc.com = lfc.com[!duplicated(lfc.com[,c('efo.id','ENTREZ')]),]

# lfc.efo = split(lfc.com, lfc.com$efo.id)
load("./data/disease.genes50.RData")
DisGen = disease.genes[same.disease == TRUE & overlap > 0 & p.adjusted < 0.05]
# remove duplicated rows, since sometimes a disease id is paired with same disease because of slight different names
DisGen = DisGen[!duplicated(DisGen$efo.id.DEGs),] #remove duplicated rows based on one column
DisGen = DisGen[,c(1,2,9)]
DisGen <- unique(DisGen %>% 
                   mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>% 
                   unnest(ensembl.id))
DisGen$ensembl.id = gsub("\"","",DisGen$ensembl.id)
DisGen$ensembl.id = gsub("c\\(","",DisGen$ensembl.id)
DisGen$ensembl.id = gsub("\\)","",DisGen$ensembl.id)
DisGen$ensembl.id = trimws(DisGen$ensembl.id)
length(unique(DisGen$efo.id))
length(unique(DisGen$ensembl.id))
DisGen$commonGenes = NULL
names(DisGen) = c("efo.id","efo.term","ensembl.id")
DisGen = merge(DisGen,lfc.com,by=c('ensembl.id','efo.id')) # merge with harmonizome

lfc.efo = split(DisGen, DisGen$efo.term)
lfc.efo = Filter(function(x) dim(x)[1] > 10, lfc.efo) # remove diseases with very few (less than 10) genes to test

rm(lfc.com)

save(lfc.efo,file="./data/disease47.genes50.lfc.RData")

##_____read log fold changes for a particular disease_______###
# load("./data/lfc.dis52.RData") 
load("./data/disease47.genes50.lfc.RData")###--get log fold change for all genes for each diseases-##

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

lfc_hgnc <- list()
for (i in 1:length(lfc.efo)) {
  lfc_hgnc[[i]] = setNames(as.numeric(lfc.efo[[i]][[4]]),as.character(lfc.efo[[i]][[6]]))
  names(lfc_hgnc)[[i]] = names(lfc.efo)[[i]]
}

##_____Create a vector with all Gene universe to proxy Array Genes_________###
hgnc_all = unique(gene.id$HGNC)
ensembl_all = unique(gene.id$ensembl.id)
entrez_all = unique(gene.id$ENTREZ)
entrezID_all = unique(gsub("^","ENTREZID:",gene.id$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_hgnc,lfc_ensembl,lfc_entrez,lfc_entrezID,hgnc_all,ensembl_all,entrez_all,entrezID_all,file="./data/disease47.genes50.lfc.namedVec.RData")
# save(lfc_hgnc,hgnc_all,file="./data/disease.genes50.lfc.namedVec.HGNC.RData")
rm(lfc.efo,gene.id)

#################################3

#####________SPIA____________#####

library(SPIA)
setwd("/home/memon/projects/msdrp/")
#___________KEGG SPIA______________#

# load("./data/lfc.namedVec.dis52.RData")
load("./data/disease47.genes50.lfc.namedVec.RData")

spia_kegg = list()
for (i in 1:length(lfc_entrez)) {
  spia_kegg[[i]] = spia(de = lfc_entrez[[i]], all = entrez_all, data.dir="./data/real_kegg/",organism="hsa")
}
names(spia_kegg) = names(lfc_entrez)
save(spia_kegg,file = "./data/spia/spia_kegg_disease47.genes50_results.RData")

# rm(list=ls())
# gc()

#load("./data/lfc.all.RData")
spia_kegg_fake = list()
for (i in 1:length(lfc_entrez)) {
  spia_kegg_fake[[i]] = spia(de = lfc_entrez[[i]], all = entrez_all, data.dir="./data/pseudo_kegg/",organism="hsa")
}
save(spia_kegg_fake,file = "./data/spia/spia_kegg_disease47.genes50_fake_results.RData")
plotP(spia_kegg[[13]])
rm(list=ls())
gc()

load("./data/spia/spia_kegg_disease.genes50_results.RData")
names(lfc_entrez$EFO_0000249)

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


dmap = read.csv("/home/memon/projects/msdrp/2725_drugs_details.csv",header = F)
dmap$V1 = gsub("[A-z]+\\':","",dmap$V1)
dmap$V1 = gsub("u","",dmap$V1)
dmap$V1 = gsub("\'","",dmap$V1)
dmap$V1 = gsub(" ","",dmap$V1)
dmap$V1 = gsub("\\{","",dmap$V1)
dmap$V1 = gsub("\\}","",dmap$V1)

library(stringr)
dmap = as.data.frame(str_split_fixed(dmap$V1, ",", 26))
dmap = dmap[,c(2,26,1,3,9)]
names(dmap) = c("chembl.id","chembl.name","phase","indication","ruleof5")
dmap= dmap %>% mutate_all(as.character)
dmap$chembl.name = gsub("[0-9]+.[0-9]+,","",dmap$chembl.name)
dmap[dmap == "None"] <- NA
save(dmap,file = "./data/drug2715details.RData")
load("./data/drug2715details.RData")


