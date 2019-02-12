#' 5th script
#' summary:
#' Download LINCS L1000 dataset 
#' Map from LINCS IDs to Chembl IDs using to PubChem IDs as intermediate 


library(EnsDb.Hsapiens.v86)
library(data.table)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

setwd("/home/memon/projects/msdrp/")

#####  get LINCS L1000 data from Harmonizome  #####

# harmonizome <- fread("gunzip -c ./data/gene_attribute_edges.txt.gz", header = TRUE, skip = 1, select = c(1, 3, 4, 7))
# http://amp.pharm.mssm.edu/Harmonizome/dataset/LINCS+L1000+CMAP+Signatures+of+Differentially+Expressed+Genes+for+Small+Molecules
# use harmonizome_downloader.ipynb for down load harmonizome dataset, OR
url ="http://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/lincscmapchemical/gene_attribute_edges.txt.gz"

if(!http_error(url) == TRUE){
  tmp <- tempfile()
  download.file(url,tmp)
  L1000 <- read.csv(gzfile(tmp),header = T, skip = 1, sep = "\t")[,c(1,3,4,7)]
  rm(tmp,url)
} else {
  print("The url is outdated, please update!")
}
L1000 <- data.table(L1000)

#####  get LINCS to PubChem mappings  #####

# lincs.mappings <- fread("../data/LINCS_LSM_Pubchem_ID_mappings.tsv") # NOT FOUND, no source mentioned
# contains "pert_id","pert_iname","pubchem_cid","SMILES" 
# lincs.mappings <- fread("./data/meta_SMILES.csv",select = c(1,3))

url = "http://maayanlab.net/SEP-L1000/downloads/meta_SMILES.csv"
if(!http_error(url) == TRUE){
  lincs.mappings <- read.csv(url)[,c(1,3)]
  rm(url)
} else {
  print("The url is outdated, please updatde!")
}
lincs.mappings <- data.table(lincs.mappings)

#####  map Entrez IDs to Ensembl  #####

L1000.genes <- L1000[, as.character(unique(GeneID))]
anno <- as.data.table(select(EnsDb.Hsapiens.v86, keys = L1000.genes, keytype = "ENTREZID", columns = c("GENEID", "ENTREZID")))
L1000 <- merge(anno[, ENTREZID := as.numeric(ENTREZID)], L1000, by.x = "ENTREZID", by.y = "GeneID", allow.cartesian = TRUE)

#####  map LINCS IDs to ChEMBL IDs  #####

unichem.url <- "https://www.ebi.ac.uk/unichem/rest/src_compound_id/"
unichem.res <- foreach(i = seq(lincs.mappings[, pubchem_cid]), .combine = rbind) %dopar% {
    lincs.id <- lincs.mappings[i, pert_id]
    pubchem.id <- lincs.mappings[i, pubchem_cid]
    chembl.id <- as.character(fromJSON(content(GET(paste0(unichem.url, lincs.mappings[i, pubchem_cid], "/22/1")), as = "text", encoding = "UTF-8")))
    if (length(chembl.id > 0) && startsWith(chembl.id, "CHEMBL")) {
        tmp <- data.table(lincs.id, pubchem.id, chembl.id)
    }
}

# annotate and tidy
L1000[, lincs.id := substr(`Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, 1, 13)]
L1000 <- merge(L1000, unichem.res, by = "lincs.id")
L1000 <- L1000[, .(ensembl.id = GENEID, gene.symbol = GeneSym, lincs.id, pubchem.id, chembl.id, perturbation = `Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit`, direction = weight)]
fwrite(L1000, "./data/L1000.tsv")
