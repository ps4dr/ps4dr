#' 5th script
#' summary:
#' 01: Download Drug Perturbed Gens Expression Profiles, LINCS L1000 dataset 
#' 02: Map from LINCS IDs to Chembl IDs using to PubChem IDs as intermediate 
#' unichem RESTful API was last accessed on 11 March, 2019. 

# source("http://bioconductor.org/biocLite.R")
# biocLite("EnsDb.Hsapiens.v86")

library(EnsDb.Hsapiens.v86)
library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)


#####################################################################
#TODO: Change to the directory where you cloned this repository
setwd("/home/memon/projects/msdrp/")
#####################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~get LINCS L1000 data from Harmonizome~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

url ="http://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/lincscmapchemical/gene_attribute_edges.txt.gz"

if(!http_error(url) == TRUE){
  tmp = tempfile()
  download.file(url,tmp)
  L1000 = read.csv(gzfile(tmp),header = T, skip = 1, sep = "\t")[,c(1,3,4,7)]
  rm(tmp,url)
} else {
  print("The url is outdated, please update!")
}
L1000 = data.table(unique(L1000))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~get LINCS to PubChem mappings~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

url = "http://maayanlab.net/SEP-L1000/downloads/meta_SMILES.csv"
if(!http_error(url) == TRUE){
  lincs.mappings = read.csv(url)[,c(1,3)]
  rm(url)
} else {
  print("The url is outdated, please updatde!")
}
lincs.mappings = data.table(unique(lincs.mappings))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~map Entrez IDs to Ensembl~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

L1000.genes = L1000[, as.character(unique(GeneID))]
anno = as.data.table(select(EnsDb.Hsapiens.v86, keys = L1000.genes, keytype = "ENTREZID", columns = c("GENEID", "ENTREZID")))
L1000 = merge(anno[, ENTREZID := as.numeric(ENTREZID)], L1000, by.x = "ENTREZID", by.y = "GeneID", allow.cartesian = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~map LINCS IDs to ChEMBL IDs~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

unichem.url = "https://www.ebi.ac.uk/unichem/rest/src_compound_id/"
unichem.res = foreach(i = seq(lincs.mappings[, pubchem_cid]), .combine = rbind) %dopar% {
    lincs.id = lincs.mappings[i, pert_id]
    pubchem.id = lincs.mappings[i, pubchem_cid]
    chembl.id = as.character(fromJSON(content(GET(paste0(unichem.url, lincs.mappings[i, pubchem_cid], "/22/1")), as = "text", encoding = "UTF-8")))
    if (length(chembl.id > 0) && startsWith(chembl.id, "CHEMBL")) {
        tmp = data.table(lincs.id, pubchem.id, chembl.id)
    }
}

L1000[, lincs.id := substr(`Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, 1, 13)]
L1000 = merge(L1000, unichem.res, by = "lincs.id")
L1000 = L1000[, .(ensembl.id = GENEID, gene.symbol = GeneSym, lincs.id, pubchem.id, chembl.id, perturbation = `Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, direction = weight)]
save(L1000, file="./data/L1000.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~L1000 Drugs~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' write all L1000 Drugs into a CSV file to use it to get their 
#' CHEMBL names and clinical tril information from CHEMBL via their API
#' using chemblid2name.ipynb script

L1000Drugs = unique(L1000[,5])
fwrite(L1000Drugs,"./data/L1000Drugs.csv",col.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


