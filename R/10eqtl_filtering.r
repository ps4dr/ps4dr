
library(data.table)
library(sqldf)
library(plyr)
library(tidyr)

setwd("/home/memon/projects/msdrp/")


stopgap = fread("./data/stopgap_snp.tsv")
length(unique(stopgap$gene.symbol))


# read efo.ids to therapeutic area mapping file
therap.area = fread("./data/opentargets.degs.therapeutic.tsv")
therap.area = unique(therap.area[,c(1,2)])

# split column with multiple elements into multiple rows 
opentargets.tas = unique(therap.area %>% 
                            mutate(therapeutic.area = strsplit(as.character(therapeutic.area), "\\|")) %>% 
                            unnest(therapeutic.area))
opentargets.tas = data.table(opentargets.tas)

relevant.tas = c("cardiovascular disease", "digestive system disease", "endocrine system disease", "hematological system disease", "immune system disease", "infectious disease", "metabolic disease", "neoplasm", "nervous system disease", "reproductive system disease", "respiratory system disease", "skeletal system disease")

neuroDis = opentargets.tas[therapeutic.area=="nervous system disease"]$efo.id

length(unique(opentargets.tas$therapeutic.area))

fisher.genes = fread("./data/fisher.genes.tsv")
fisher.genes = merge(fisher.genes, opentargets.tas[, .(efo.id, therapeutic.area)], by.x = "efo.id.stopgap", by.y = "efo.id", all = FALSE)
fisher.genes = unique(fisher.genes[,-4]) # give false counts if this column remains
fisher.genes = unique(fisher.genes)
dis2therap = fisher.genes[same.disease == TRUE & therapeutic.area %in% relevant.tas, .(number = sum(p.adjusted < 0.05)), by = therapeutic.area]

library(plyr)
dis2therap = count(fisher.genes,c("therapeutic.area","efo.id.stopgap"))

stopgap = unique(fread("./data/stopgap_snp.tsv"))
dis_genes = sqldf('SELECT "efo.id", count("gene.symbol") AS genes FROM stopgap GROUP BY "efo.id"')


stopgap.ad = stopgap[efo.id=="EFO_0000249"]
stopgap.pd = stopgap[efo.id=="EFO_0002508"]
save(stopgap.ad,stopgap.pd,file="./data/testadpd.RData")

stopgap.ad$snp.ld=trimws(stopgap.ad$snp.ld)             
             
load("/home/memon/genetic_analyses/scripts/pathburden/snp_gene_eqtl.RData")
mysnp = snp_gene[,1]

# xm = gwas.ld.dhscor
# xm =var2gene.vepSimp
# rm(gwas.ld.dhscor)
# xm$Gene = gsub("^$", NA, xm$Gene)
# xm = xm[which(!is.na(xm$Gene)),]
# 
# length(unique(xm$SNP.ld))
# length(unique(stopgap$snp.ld))
# length(intersect(xm$SNP.ld,stopgap$snp.ld))
# 
# length(intersect(snp_gene$RefSNP_id,stopgap$snp.ld))
# 
# length(unique(stopgap.ad$snp.ld))
# length(intersect(snp_gene$RefSNP_id,stopgap.ad$snp.ld))
# 
# length(unique(stopgap.pd$snp.ld))
# length(intersect(snp_gene$RefSNP_id,stopgap.pd$snp.ld))

stopgap.neuro = list()
k=0
for(i in neuroDis){ 
  k = k+1
  stopgap.neuro[[i]] = subset(stopgap, efo.id==i)
}
stopgap.neuro = stopgap.neuro[sapply(stopgap.neuro, function(x) dim(x)[1]) > 0]
stopgap.neuro = do.call(rbind, stopgap.neuro) #list to df

neuro = merge(stopgap.neuro,snp_gene,by.x="snp.ld",by.y="RefSNP_id")

stopgap.ad = neuro[efo.id=="EFO_0000249"]
stopgap.pd = neuro[efo.id=="EFO_0002508"]

neuro = list()
k=0
for (i in stopgap.neuro){
  k = k+1
  neuro[[i]] = subset(stopgap.neuro, snp.ld==mysnp)
}

intersect(stopgap.neuro[[1]][[1]],mysnp)

intersect(stopgap.ad$snp.ld,stopgap.pd$snp.ld)

stopgap.bestld = data.table(stopgap.bestld)
test = stopgap.bestld[snp.ld=="rs169201"]


mysnp = scan("/home/memon/projects/drug_target_prioritization/STOPGAP/stopgap2_SNPs.txt",character())
mysnp <- mysnp[!is.na(mysnp)]
mysnp <- mysnp[!is.null(mysnp)]
mysnp <- trimws(mysnp,which = "both")

library(haploR) 
mysnp3 = mysnp[100000:150000]

haplo_res_big <- as.data.frame(queryHaploreg(query = mysnp3,ldThresh = 0.8))


# get Drugbank2Dataframe
# source("/home/memon/projects/drug_target_prioritization/GCMap/dbxml_handle.R")
# drugbank = dbxml2df("/home/memon/projects/drug_target_prioritization/GCMap/data/full database.xml")
# save(drugbank,file = "/home/memon/projects/drug_target_prioritization/GCMap/data/drugbank.RData")
# write.csv(drugbank,file = "/home/memon/projects/drug_target_prioritization/GCMap/data/drugbank.csv",quote = F)
drugbank[1000,]$`external-links`

drugbank.chembl = subset(drugbank, grepl("ChEMBL",`external-identifiers`))
drugbank.chembl$`external-identifiers` = gsub(".*ChEMBL","",drugbank.chembl$`external-identifiers`)
