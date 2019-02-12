library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(doMC)
no_cores <- detectCores() - 4

library(data.table) # (V 1.11.6)
library(dplyr) # (V 0.7.6)
library(tidyr) # (V 0.8.1)
library(SPIA)

setwd("/home/memon/projects/msdrp/")

harmonizome <- unique(fread("./dat/harmonizome.tsv"))
harmonizome$perturbation <- NULL
harmonizome = merge(harmonizome,gene.id.entrez,by="ensembl.id")
harmonizome = unique(harmonizome[,c(1,7,5,6)])
harmonizome = harmonizome[!duplicated(harmonizome[,c('chembl.id', 'ensembl.id')]),]


load("./dat/super.drugs.48D.padj1e-5.RData")
super.drugs = super.drugs[,c(1,2,8)]

drug.genes <- unique(super.drugs %>% 
                       mutate(ensembl.id = strsplit(as.character(commonGenes), ",")) %>% 
                       unnest(ensembl.id))
drug.genes$ensembl.id = gsub("\"","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("c\\(","",drug.genes$ensembl.id)
drug.genes$ensembl.id = gsub("\\)","",drug.genes$ensembl.id)
drug.genes$ensembl.id = trimws(drug.genes$ensembl.id)

length(unique(drug.genes$chembl.id)) # 2290
length(unique(drug.genes$ensembl.id)) # 1790

drug.genes = merge(drug.genes,harmonizome,by=c('chembl.id','ensembl.id')) # merge with harmonizome
drug.genes = drug.genes[!duplicated(drug.genes[,c('chembl.id','ensembl.id','efo.id')]),]
length(unique(drug.genes$efo.id))

fd.disease = split(drug.genes, drug.genes$efo.id) # split on diseases first
fd.disease.drug = lapply(fd.disease, function(x) split(x, x$chembl.id)) # split on drug inside each disease
#m = fd.disease.drug[lapply(fd.disease.drug, length) > 2]

### remove chemble id and efo.id columns sicne they are not required any more
for (element in 1:length(fd.disease.drug)){
  for (subelem in 1:length(fd.disease.drug[[element]])){
    fd.disease.drug[[element]][[subelem]]$chembl.id = NULL
    fd.disease.drug[[element]][[subelem]]$ensembl.id = NULL
    fd.disease.drug[[element]][[subelem]]$efo.id = NULL
  }
}
save(fd.disease.drug,file="./dat/fd.disease38.drug.RData")

##_____Create Named Vector for log fold changes in each disease_____________###

load("./dat/fd.disease38.drug.RData")
load("./dat/gene.universe.RData")

lfc_drug_entrez = fd.disease.drug
for (element in 1:length(fd.disease.drug)){
  for (subelem in 1:length(fd.disease.drug[[element]])){
    lfc_drug_entrez[[element]][[subelem]] = setNames(as.numeric(fd.disease.drug[[element]][[subelem]][[2]]),as.character(fd.disease.drug[[element]][[subelem]][[1]]))
  }
}

lfc_drug_entrezID = fd.disease.drug
for (element in 1:length(fd.disease.drug)){
  for (subelem in 1:length(fd.disease.drug[[element]])){
    lfc_drug_entrezID[[element]][[subelem]] = setNames(as.numeric(fd.disease.drug[[element]][[subelem]][[2]]),as.character(gsub("^","ENTREZID:",fd.disease.drug[[element]][[subelem]][[1]])))
  }
}
save(lfc_drug_entrez,lfc_drug_entrezID,entrez_all,entrezID_all,file="./dat/lfc.dis38.drug_gene.big.RData")

##---SET working Directory---###
setwd("/home/memon/projects/drugrep/drug2path")

##---Main Function---###
spia_fun <- function(x){
  spia_drug_kegg = list()
  spia_drug_kegg = spia(de = x, all = entrez_all, data.dir="./dat/real_kegg/",organism="hsa")
}

##---Get Disease-Drugs list---###
load("./dat/lfc.dis38.drug_gene.big.RData")
lfc_test = lfc_drug_entrez[3:4]
rm(lfc_drug_entrez,lfc_drug_entrezID,entrezID_all)

##____Single Node______###

spia_kegg_2 = list()

for (diseases in 1:length(lfc_test)){
  spia_kegg_2[[diseases]] <- mclapply(lfc_test[[diseases]], spia_fun, mc.cores = no_cores)
  names(spia_kegg_2)[[diseases]] = names(lfc_test)[[diseases]]
}

save(spia_kegg_2,file = "./dat/spia/spia_kegg_38_nopar_test1.RData")


