#' 6th script
#' summary:

library(stringr)
library(SPIA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####_____DEGs Data Preperation for SPIA______#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

# load("./data/drugPdisease.genes.RData")
# drugPdisease.genes = drugPdisease.genes[,c(1,2,9)]
load("./data/drugPdisease.genes.48D.padj.RData")
load("./data/drugPdisease.genes.48D.padj1e-5.RData")
drugPdisease.genes = drugPdisease.genes[,c(1,2,8)]

drugGWAS.genes = unique(fread("./data/fisher.drugs.commongenes.tsv"))
drugGWAS.genes = drugGWAS.genes[,c(1,3,9)]
tmp =  drugGWAS.genes
#drugPdisease.genes.chembl = split(drugPdisease.genes, drugPdisease.genes$chembl.id)

drug.genes <- unique(tmp %>% 
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



disease.genes <- fread("./data/disease.genes.commongenes.tsv")
DisGen = disease.genes[same.disease == TRUE & overlap > 0]
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


dmap = read.csv("/home/memon/projects/msdrp/2725_drugs_details.csv",header = F)
dmap$V1 = gsub("[A-z]+\\':","",dmap$V1)
dmap$V1 = gsub("u","",dmap$V1)
dmap$V1 = gsub("\'","",dmap$V1)
dmap$V1 = gsub(" ","",dmap$V1)
dmap$V1 = gsub("\\{","",dmap$V1)
dmap$V1 = gsub("\\}","",dmap$V1)

dmap = as.data.frame(str_split_fixed(dmap$V1, ",", 26))
dmap = dmap[,c(2,26,1,3,9)]
names(dmap) = c("chembl.id","chembl.name","phase","indication","ruleof5")
dmap= dmap %>% mutate_all(as.character)
dmap$chembl.name = gsub("[0-9]+.[0-9]+,","",dmap$chembl.name)
dmap[dmap == "None"] <- NA
save(dmap,file = "./data/drug2715details.RData")
load("./data/drug2715details.RData")


