#' 3rd script
#' summary:
#' Mapping of MeSH terms to EFO IDs from GWASs Data processed by STOPGAP Pipeline
#' for downloading Differential Gene Expression Data from Open Targtets by EFO IDs
#' EFO ontology (version 2.105) file has been used for these mapping purpose

library(EnsDb.Hsapiens.v86)
library(data.table)
library(dplyr)
library(tidyr)
library(httr)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
resultsFolder = args[1]
if (! file.exists(resultsFolder)) {
    dir.create(resultsFolder)
}
print("Using results folder at " + resultsFolder)

# get GWASs Data from STOPGAP pipeline output 
load("./data/stopgap.gene.mesh.RData")
GWASs = stopgap.gene.mesh
rm(stopgap.gene.mesh)

length(unique(GWASs$msh))
length(unique(GWASs$gene.v19))
length(unique(GWASs$snp.ld))

# keep relevant columns
GWASs = data.table(GWASs)
GWASs = unique(GWASs[, .(snp.ld,gene.v19, msh, pvalue, gene.score, gene.rank.min, source)])
# replace p-values of zero with arbitrarily low p-value
GWASs[pvalue == 0, pvalue := 3e-324]

# map gene symbols to Ensembl
GWASs.genes = unique(GWASs[, gene.v19])
anno = as.data.table(select(EnsDb.Hsapiens.v86, keys = GWASs.genes, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL")))
GWASs = merge(GWASs, anno, by.x = "gene.v19", by.y = "SYMBOL", all = FALSE)

#' read EFO ontolgy (version 2.105) files in csv format downloaded from bioportal on 11th March, 2019
#' "http://data.bioontology.org/ontologies/EFO/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv"

efo = fread("/home/memon/projects/msdrp/EFO.csv")
efo = efo[,c(1,2,3)] # filter out unnecessary columns
efo$`Class ID` = gsub(".*efo/","",efo$`Class ID`) # replace urls to get only EFO IDs
efo = efo[efo$`Class ID` %like% "EFO",] # filter out other disease IDs other than EFO
efo = unite(efo, Disease, c('Preferred Label', Synonyms), sep = "|",remove=FALSE)
efo = efo[,c(1,2)]

# create new rows for each synonyms of EFO IDs 
efo = unique(efo %>% 
                mutate(Disease = strsplit(as.character(Disease), "\\|")) %>% 
                unnest(Disease))
efo$Disease = gsub("^$", NA, efo$Disease)
efo = efo[which(!is.na(efo$Disease)),]
efo$upper = toupper(efo$Disease)
efo$upper = trimws(efo$upper, which = "both")

mesh.terms = unique(GWASs[,3])
mesh.terms$upper = toupper(mesh.terms$msh)
mesh.terms$upper = trimws(mesh.terms$upper, which = "both")

# merge MeSH terms from GWASs to efo IDs
mesh2efo = merge(mesh.terms,efo, by="upper")
mesh2efo = mesh2efo[,c(2,3,4)]
names(mesh2efo) = c("mesh","efo.id","efo.term")

# merge GWASs with mesh2efo table
GWASs = merge(GWASs, mesh2efo, by.x = "msh", by.y = "mesh", all = FALSE)
GWASs = GWASs[, .(snp.ld,ensembl.id = GENEID, gene.symbol = gene.v19, efo.id, efo.term, GWASs.pvalue = pvalue, GWASs.gene.score = gene.score, GWASs.gene.rank = gene.rank.min,source)]
unique(GWASs$efo.id)

save(GWASs, file="./data/GWASs.RData")
length(unique(GWASs$efo.id))
length(unique(GWASs$ensembl.id))
stopgap = fread("./data/stopgap.tsv")
length(unique(stopgap$efo.id))
length(unique(stopgap$ensembl.id))
