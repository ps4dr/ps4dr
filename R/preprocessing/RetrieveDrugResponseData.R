#' 3rd script
#' summary:
#' 01: Download Drug Perturbed Gens Expression Profiles, LINCS L1000 dataset 
#' 02: Map from LINCS IDs to Chembl IDs using to PubChem IDs as intermediate 
#' unichem RESTful API was last accessed on 11 March, 2019. 


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(httr)))
suppressWarnings(suppressMessages(library(jsonlite)))

#####################################################################
#TODO: Change to the directory where you cloned this repository
#~~~~~~~Using relative path~~~~~~~#
ensureFolder = function(folder) {
  if (! file.exists(folder)) {
    dir.create(folder)
  }
}

args = commandArgs(trailingOnly = TRUE)
resultsFolder = normalizePath(args[1])
ensureFolder(resultsFolder)
sprintf("Using results folder at %s", resultsFolder)

dataFolder = file.path(resultsFolder)
#####################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~get LINCS L1000 data from Harmonizome~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!file.exists(file.path(dataFolder, "L1000_raw.RData"))){
  url ="http://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/lincscmapchemical/gene_attribute_edges.txt.gz"
  tryCatch(if(!http_error(url) == TRUE){
    tmp = tempfile()
    download.file(url,tmp)
    L1000_raw = read.csv(gzfile(tmp),header = T, skip = 1, sep = "\t")[,c(1,3,4,7)]
    rm(tmp,url)
  } else {
    print("The url is outdated, please update!")
  },
  error=function(e) 1)
  L1000_raw = data.table(unique(L1000_raw))
  L1000_raw[, lincs_id := substr(`Perturbation.ID_Perturbagen_Cell.Line_Time_Time.Unit_Dose_Dose.Unit`, 1, 13)]
  L1000_raw$Perturbation.ID_Perturbagen_Cell.Line_Time_Time.Unit_Dose_Dose.Unit = NULL
  L1000_raw = data.table(unique(L1000_raw))
  save(L1000_raw,file = file.path(dataFolder, "L1000_raw.RData"))
} else {cat(sprintf("~~ L1000_raw file already exists, not downloading again. ~~\n"))
  load(file.path(dataFolder, "L1000_raw.RData"))}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~get LINCS to PubChem mappings~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!file.exists(file.path(dataFolder, "lincs_pubchem.RData"))){
  url ="http://maayanlab.net/SEP-L1000/downloads/meta_SMILES.csv"
  tryCatch(if(!http_error(url) == TRUE){
    cat("~~Downloading lincs_pubchem mapping file.~~")
    lincs_pubchem = read.csv(url)[,c(1,3)]
    rm(url)
  } else {
    cat("The url is outdated, please updatde!")
  },
  error=function(e) 1)
  
  lincs_pubchem = lincs_pubchem[which(! is.na(lincs_pubchem$pubchem_cid)),]
  lincs_pubchem = data.table(unique(lincs_pubchem))
  lincs_pubchem$pert_id = as.character(lincs_pubchem$pert_id)
  lincs_pubchem$pubchem_cid = as.character(lincs_pubchem$pubchem_cid)
  save(lincs_pubchem,file = file.path(dataFolder, "lincs_pubchem.RData"))
} else {cat(sprintf("~~ lincs_pubchem file already exists, not downloading again. ~~\n"))
  load(file.path(dataFolder, "lincs_pubchem.RData"))}

lincs_pubchem = merge(L1000_raw, lincs_pubchem, by.x= "lincs_id",by.y = "pert_id")
lincs_pubchem = unique(lincs_pubchem[,.(lincs_id,pubchem_cid)])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~map Entrez IDs to Ensembl~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
load(file.path(dataFolder, "geneID_97.RData"))

L1000_genes = L1000_raw[, as.character(unique(GeneID))]
L1000 = merge(L1000_raw, gene_id, by.x = "GeneID", by.y = "ENTREZ")
L1000 = L1000[,c(1,5,2,3,4)]
names(L1000) = c("ENTREZID","GENEID","GeneSym","weight","lincs_id")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~map LINCS IDs to ChEMBL IDs~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!file.exists(file.path(dataFolder, "L1000_v97.RData"))){
  cat(sprintf("~~ Mapping BROAD IDs to ChEMBL via PubChem IDs. ~~\n"))
  pb <- txtProgressBar(min = 0, max = length(lincs_pubchem[, pubchem_cid]), style = 3)
  unichem_url = "https://www.ebi.ac.uk/unichem/rest/src_compound_id/"
  unichem_map = data.table()
  tryCatch(for(i in 1:length(lincs_pubchem[, pubchem_cid])){
    Sys.sleep(0.1)
    lincs_id = lincs_pubchem[i, lincs_id]
    pubchem_id = lincs_pubchem[i, pubchem_cid]
    chembl_id = as.character(fromJSON(content(GET(paste0(unichem_url, lincs_pubchem[i, pubchem_cid], "/22/1")), as = "text", encoding = "UTF-8")))
    if (length(chembl_id > 0) && startsWith(chembl_id, "CHEMBL")) {
      tmp = data.table(lincs_id, pubchem_id, chembl_id)
      unichem_map = rbind(unichem_map,tmp)
    }
    setTxtProgressBar(pb, i)
  },
  error=function(e) 1)
  close(pb)
  
  L1000 = merge(L1000, unichem_map, by = "lincs_id")
  L1000 = L1000[, .(ensembl.id = GENEID, gene.symbol = GeneSym, lincs.id=lincs_id, pubchem.id=pubchem_id, chembl.id=chembl_id, direction = weight)]
  
  save(L1000, file=file.path(dataFolder, "L1000_v97.RData"))
  
} else {cat(sprintf("~~ L1000 file already exists, not mapping again. ~~\n"))
  load(file.path(dataFolder, "L1000_v97.RData"))}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~L1000 Drugs~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' write all L1000 Drugs into a CSV file to use it to get their 
#' CHEMBL names and clinical tril information from CHEMBL via their API
#' using chemblid2name.ipynb script

L1000Drugs = unique(L1000[, 5])
fwrite(L1000Drugs, file=file.path(dataFolder,"L1000Drugs.csv"), col.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

