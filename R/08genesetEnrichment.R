#' 1st script
#' summary:
#' cleaning up STOPGAP data


#source("https://bioconductor.org/biocLite.R")
#biocLite("EnrichmentBrowser")
#biocLite("SPIA")
#biocLite("graphite")

library(data.table)
library(tidyverse)
library(SPIA)
library(graphite)

setwd("/home/memon/projects/msdrp/")


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

#__map ensembl IDs to ENTREZ ID___#

#-----get mapping among ENTREZ_HGNC_ENSEMBL_IDs---------#
# library(biomaRt)
# ensembl = useEnsembl(biomart="ensembl", version=92, dataset="hsapiens_gene_ensembl") # v92 has less id than v79
# val <- c(1:23,"X","Y")
# gene.id <- getBM(attributes=c('entrezgene','hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position'),
#                    filters ='chromosome_name', values =val, mart = ensembl)
# 
# save(gene.id,file = "./dat/geneID.RData")
# rm(ensembl,val)

# load("./dat/geneID.RData")
# gene.id$entrezgene <- gsub("^$", NA, gene.id$entrezgene)
# gene.id <- gene.id[which(!is.na(gene.id$entrezgene)),]
# gene.id <- data.table(gene.id)
# gene.id<- unique(gene.id[,c('entrezgene','ensembl_gene_id')])
# names(gene.id) <- c("ENTREZ","ensembl.id")
# gene.id <- gene.id[!duplicated(gene.id$ensembl.id),]
# gene.id.entrez <- gene.id[!duplicated(gene.id$ensembl.id),]
# save(gene.id.entrez,file = "./dat/gene.id.entrez.RData")

load("./dat/gene.id.entrez.RData")
lfc.com <- merge(lfc.com,gene.id.entrez,by="ensembl.id")
lfc.com = lfc.com[!duplicated(lfc.com[,c('efo.id','ENTREZ')]),]

lfc.efo = split(lfc.com, lfc.com$efo.id)
rm(degs,lfc.neg,lfc.pos,lfc.com)

#degs.ad <- lfc.com[efo.id == "EFO_0000249"]

#save(lfc.com,file="./dat/lfc.com.RData")
save(lfc.efo,file="./dat/lfc.efo.RData")

##_____read log fold changes for a particular disease_______###
load("./dat/lfc.efo.RData") ###--get log fold change for all genes for each diseases-##

##_____Create Named Vector for log fold changes in each disease_____________###

lfc_ensembl <- list()
for (i in 1:length(lfc.efo)) {
  lfc_ensembl[[i]] = setNames(as.numeric(lfc.efo[[i]][[3]]),as.character(lfc.efo[[i]][[1]]))
}

lfc_entrez <- list()
for (i in 1:length(lfc.efo)) {
  lfc_entrez[[i]] = setNames(as.numeric(lfc.efo[[i]][[3]]),as.character(lfc.efo[[i]][[4]]))
}

lfc_entrezID <- list()
for (i in 1:length(lfc.efo)) {
  lfc_entrezID[[i]] = setNames(as.numeric(lfc.efo[[i]][[3]]),as.character(gsub("^","ENTREZID:",lfc.efo[[i]][[4]])))
}


##_____Create a vector with all Gene universe to proxy Array Genes_________###
ensembl_all = unique(gene.id.entrez$ensembl.id)
entrez_all = unique(gene.id.entrez$ENTREZ)
entrezID_all = unique(gsub("^","ENTREZID:",gene.id.entrez$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_ensembl,lfc_entrez,lfc_entrezID,ensembl_all,entrez_all,entrezID_all,file="./dat/lfc.all.RData")
rm(lfc.efo,gene.id.entrez)

#################################3

#####__Preapre SPIA Pathway Data set__####

##_____KEGG SPIA____###
##----------Prepare SPIA from KEGG pathway kgml files---------##
# makeSPIAdata(kgml.path="./dat/kgml/",organism="hsa",out.path="./dat/real_kegg/")

##_____REACT SPIA____###

# reactome <- pathways("hsapiens", "reactome")
# reactome <- convertIdentifiers(reactome, "ENTREZID")
# save(reactome,file = "./dat/reactomedb.RData")
# load("./dat/reactomedb.RData") # better load it here, sicne it takes long to map reactome default uniprot id to ENTREZ ID
# prepareSPIA(reactome, "./dat/real_react/hsa")

#####______BioCarta____________#####
# biocarta <- pathways("hsapiens", "biocarta")
# biocarta <- convertIdentifiers(biocarta, "ENTREZID")
# prepareSPIA(biocarta, "./dat/real_biocarta/hsa")

##_______________________________###


###___Create Random (Pseudo) Pathways____####

# #___create universe with ENTREZID:____#
# load("./dat/geneID.RData")
# universe = unique(na.omit(gene.id.entrez$entrezgene)) # for kegg
# universe = unique(na.omit(gsub("^","ENTREZID:",gene.id.entrez$entrezgene))) # for reactome, biocarta
# 
# load("./dat/real_kegg/hsaSPIA.RData")
# load("./dat/real_react/hsaSPIA.RData")
# load("./dat/real_biocarta/hsaSPIA.RData")
# 
# pseudo_path = path.info
# rm(path.info)
# 
# #_____Shuffle first 25 matrices in pseudo_path______#
# for (mylist in 1:length(pseudo_path)) {
#   for (mymat in 1:(length(pseudo_path[[mylist]])-3)) { #not looping for last 3 indices since they are not matrices
#     pseudo_path[[mylist]][[mymat]][] = pseudo_path[[mylist]][[mymat]][sample(length(pseudo_path[[mylist]][[mymat]]))]
#   }
# }
# 
# #___Replace ENTREZ ID with ENTREZ ID from Universe_____#
# for (mylist in 1:length(pseudo_path)) {
#   for (mymat in length(pseudo_path[[mylist]])-1) { #need to change index number to '1' for reactome and biocarta, and to '2' for kegg
#     pseudo_path[[mylist]][[mymat]] = as.character(sample(universe,size = length(pseudo_path[[mylist]][[mymat]]),replace = TRUE))
#   }
# }
# 
# #_____change row & column names of (25) matrices______#
# for (mylist in 1:length(pseudo_path)) {
#   for (mymat in 1:(length(pseudo_path[[mylist]])-3)) { #not looping for last 3 indices since they are not matrices
#     colnames(pseudo_path[[mylist]][[mymat]]) = pseudo_path[[mylist]][[27]] #change index to '27' for reactome and biocarta, and to '26' for kegg
#     rownames(pseudo_path[[mylist]][[mymat]]) = pseudo_path[[mylist]][[27]] #change index to '27' for reactome and biocarta, and to '26' for kegg
#   }
# }
# 
# path.info <- pseudo_path
# 
# save(path.info,file = "./dat/pseudo_kegg/hsaSPIA.RData")
# save(path.info,file = "./dat/pseudo_react/hsaSPIA.RData")
# save(path.info,file = "./dat/pseudo_biocarta/hsaSPIA.RData")
# rm(path.info,pseudo_path)
# #_______________________________________________________##





# #__________________________#
#####___Drug Data Preperation___#####

load("./dat/gene.id.entrez.RData")
fisher.drugs = fread("./dat/fisher.drugs.commongenes.tsv")
fisher.drugs = as.data.frame(unique(fisher.drugs$chembl.id))
names(fisher.drugs) = "chembl.id"

harmonizome = unique(fread("./dat/harmonizome.tsv"))
harmonizome = harmonizome[,c(5,1,7)]
harmonizome = harmonizome[order(ensembl.id,decreasing = TRUE), ]
harmonizome = harmonizome[!duplicated(harmonizome[,c('ensembl.id', 'chembl.id')]),]
harmonizome = merge(harmonizome,gene.id.entrez,by="ensembl.id")
harmonizome = merge(fisher.drugs,harmonizome,by="chembl.id")
harmonizome = harmonizome[,c(1,2,4,3)]
#harmonizome.chembl = split(harmonizome, harmonizome$chembl.id)

load("./dat/super.drugs.RData")
super.drugs = super.drugs[,c(1,2,9)]
super.drugs.chembl = split(super.drugs, super.drugs$chembl.id)

drug.genes <- unique(super.drugs %>% 
                       mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>% 
                       unnest(ensembl.id))
drug.genes$ensembl.id = gsub("\"","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("c\\(","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("\\)","",drug.genes$ensembl.id)
drug.genes$ensembl.id = trimws(drug.genes$ensembl.id)
drug.genes = merge(drug.genes,harmonizome,by=c('chembl.id','ensembl.id')) # merge with harmonizome

drug.genes.list <- split(drug.genes, list(drug.genes$chembl.id,drug.genes$efo.id)) #
drug.genes.list= Filter(function(x) dim(x)[1] > 0, drug.genes.list) # remove empty lists
save(drug.genes.list,file = "./dat/drug.genes.list.RData")

# length(intersect(harmonizome.chembl[["CHEMBL98"]][["ensembl.id"]],drug.genes.list[["CHEMBL98.EFO_0000249"]][["ensembl.id"]]))
# length(harmonizome.chembl[["CHEMBL98"]][["ensembl.id"]])
# length(drug.genes.list[["CHEMBL98.EFO_0000249"]][["ensembl.id"]])

####____Create named entity list for each drug_disease pairs___####
load("./dat/drug.genes.list.RData")
lfc_drug_ensembl <- list()
for (i in 1:length(drug.genes.list)) {
  lfc_drug_ensembl[[i]] = setNames(as.numeric(drug.genes.list[[i]][[5]]),as.character(drug.genes.list[[i]][[2]]))
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
ensembl_all = unique(gene.id.entrez$ensembl.id)
entrez_all = unique(gene.id.entrez$ENTREZ)
entrezID_all = unique(gsub("^","ENTREZID:",gene.id.entrez$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_drug_ensembl,lfc_drug_entrez,lfc_drug_entrezID,ensembl_all,entrez_all,entrezID_all,file="./dat/lfc.drug.all.RData")

##____________________________________________________##


