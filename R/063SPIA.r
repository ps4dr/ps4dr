library(data.table)
library(dplyr)
library(tidyr)

setwd("/home/memon/projects/msdrp/")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####___Drug Data Preperation for SPIA___#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

load("./data/drug2715details.RData")
dmap = data.table(dmap[,c(1,2,3)])
dmap = dmap[phase==4|phase==3|phase==2|phase==1]
harmonizome = merge(dmap[,c(1,2)],harmonizome,by="chembl.id")

# load("./data/super.drugs.RData")
# super.drugs = super.drugs[,c(1,2,9)]
# load("./data/super.drugs.48D.padj1e-5.RData")

##----Drug Perturbed Disease Genes (DEGs_GWAS)---##
# load("./data/drugPdisease.genes50.padj1e-5.RData")
# super.drugs = drugPdisease.genes[,c(2,4,9)]

##----Drug Perturbed GWAS Genes---##
load("./data/drugGWAS.genes50.padj1e-10.RData")
super.drugs = drugGWAS.genes[,c(2,4,9)]

#super.drugs.chembl = split(super.drugs, super.drugs$chembl.id)

drug.genes <- unique(super.drugs %>% 
                       mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>% 
                       unnest(ensembl.id))
drug.genes$ensembl.id = gsub("\"","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("c\\(","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("\\)","",drug.genes$ensembl.id)
drug.genes$ensembl.id = trimws(drug.genes$ensembl.id)

length(unique(drug.genes$chembl.name))
length(unique(drug.genes$efo.term))
length(unique(drug.genes$ensembl.id))

drug.genes = merge(drug.genes,harmonizome,by=c('chembl.name','ensembl.id')) # merge with harmonizome
drug.genes = drug.genes[!duplicated(drug.genes[,c('chembl.id','ensembl.id','efo.term')]),]

# drug.genes.list <- split(drug.genes, list(drug.genes$chembl.id,drug.genes$efo.id)) #
# drug.genes.list <- split(drug.genes, list(drug.genes$chembl.name,drug.genes$efo.term))
# fd.disease.drug = Filter(function(x) dim(x)[1] > 0, fd.disease.drug) # remove empty lists

fd.disease = split(drug.genes, drug.genes$efo.term) # split on diseases first
fd.disease.drug = lapply(fd.disease, function(x) split(x, x$chembl.name)) # split on drug inside each disease
fd.disease.drug = fd.disease.drug[lapply(fd.disease.drug, length) > 10] # remove diseases with very few drugs to test


### remove chemble.name, chemble.id and efo.term columns sicne they are not required any more
for (element in 1:length(fd.disease.drug)){
  for (subelem in 1:length(fd.disease.drug[[element]])){
    fd.disease.drug[[element]][[subelem]]$chembl.name = NULL
    fd.disease.drug[[element]][[subelem]]$chembl.id = NULL
    fd.disease.drug[[element]][[subelem]]$efo.term = NULL
  }
}
# save(drug.genes.list,file = "./data/drug.genes.list.48D.padj.RData")
# save(fd.disease.drug,file = "./data/drug.genes.list.drugPdisease.genes50.padj1e-5.RData")
save(fd.disease.drug,file = "./data/drug.genes.list.drugGWAS.genes50.padj1e-10.RData")

####____Create named entity list for each drug_disease pairs___####
# load("./data/drug.genes.list.RData")
# load("./data/drug.genes.list.drugPdisease.genes50.padj1e-5.RData")
load("./data/drug.genes.list.drugGWAS.genes50.padj1e-10.RData")

lfc_drug_ensembl <- fd.disease.drug
for (i in 1:length(fd.disease.drug)) {
  for (j in 1:length(fd.disease.drug[[i]])){
    lfc_drug_ensembl[[i]][[j]] = setNames(as.numeric(fd.disease.drug[[i]][[j]][[2]]),as.character(fd.disease.drug[[i]][[j]][[1]]))
  }
}

lfc_drug_entrez <- fd.disease.drug
for (i in 1:length(fd.disease.drug)) {
  for (j in 1:length(fd.disease.drug[[i]])){
    lfc_drug_entrez[[i]][[j]] = setNames(as.numeric(fd.disease.drug[[i]][[j]][[2]]),as.character(fd.disease.drug[[i]][[j]][[3]]))
    #names(lfc_drug_entrez)[i] = names(fd.disease.drug)[i]
    #names(lfc_drug_entrez[[i]][j]) = names(fd.disease.drug[[i]][j])
  }
}

lfc_drug_entrezID <- fd.disease.drug
for (i in 1:length(fd.disease.drug)) {
  for (j in 1:length(fd.disease.drug[[i]])){
    lfc_drug_entrezID[[i]][[j]] = setNames(as.numeric(fd.disease.drug[[i]][[j]][[2]]),as.character(gsub("^","ENTREZID:",fd.disease.drug[[i]][[j]][[3]])))
  }
}

##_____Create a vector with all Gene universe to proxy Array Genes_________###
load("./data/gene.id.entrez.RData")
ensembl_all = unique(gene.id.entrez$ensembl.id)
entrez_all = unique(gene.id.entrez$ENTREZ)
entrezID_all = unique(gsub("^","ENTREZID:",gene.id.entrez$ENTREZ)) ## with ENTREZID: in front of each id

save(lfc_drug_ensembl,ensembl_all,file="./data/lfc.drug.drugGWAS.genes50.padj1e-10.ensembl.RData")
save(lfc_drug_entrez,entrez_all,file="./data/lfc.drug.drugGWAS.genes50.padj1e-10.entrez.RData")
save(lfc_drug_entrezID,entrezID_all,file="./data/lfc.drug.drugGWAS.genes50.padj1e-10.entrezID.RData")

# save(lfc_drug_ensembl,ensembl_all,file="./data/lfc.drug.drugPdisease.genes50.padj1e-5.ensembl.RData")
# save(lfc_drug_entrez,entrez_all,file="./data/lfc.drug.drugPdisease.genes50.padj1e-5.entrez.RData")
# save(lfc_drug_entrezID,entrezID_all,file="./data/lfc.drug.drugPdisease.genes50.padj1e-5.entrezID.RData")


##____________________________________________________##


