#' This script creates pathway files for SPIA
#' from excel files

###~~~~~main function to read all sheets in a excel file~~~~~###
read_allsheets <- function(filename, tibble = FALSE) {
  require(readxl)
  sheets <- readxl::excel_sheets(filename)
  tmp <- lapply(sheets, function(x) readxl::read_excel(filename, sheet = x,col_types ="numeric",col_names = TRUE))
  if(!tibble) tmp <- lapply(tmp, as.matrix)
  tmp <- lapply(tmp,function(x) {rownames(x) <- colnames(x); x})
  names(tmp) <- sheets
  tmp
}

##~~~read all excel fles in a folder to create te main lsit objects~~~~~##
##~~~keep all excel files in a single folder~~~~~~##
##~~~ change the directory accordingly each time 

# setwd("/home/memon/projects/msdrp/kegg/")
# setwd("/home/memon/projects/msdrp/reactome/")
setwd("/home/memon/projects/msdrp/wikipath/")

filenames <-list.files(path = ".")
path.info = lapply(paste0(getwd(),"/", filenames), read_allsheets) 
filenames = gsub('_unflatten.xlsx','',filenames)
filenames = gsub('hsa','',filenames)
filenames = gsub('R-HSA-','',filenames)
filenames = gsub('.xlsx','',filenames)
names(path.info) = filenames

###~~~~~~add last 3 objects in each list inside the path.info~~~~###

for (pathway in 1:length(path.info)) {
  for (matrix in 1:length(path.info[pathway])) {
    path.info[[pathway]][['nodes']] <- rownames(path.info[[pathway]][[matrix]])
    path.info[[pathway]][['title']] <- filenames[pathway]
    path.info[[pathway]][['NumberOfReactions']] <- as.integer(0)
    names(path.info[[pathway]])[[3]] <- "binding/association"
    names(path.info[[pathway]])[[22]] <- "activation_binding/association"
    
  }
}

##############################################
# save(path.info,file="/home/memon/projects/msdrp/kegg/hsaSPIA.RData")
# save(path.info,file="/home/memon/projects/msdrp/reactome/hsaSPIA.RData")
# save(path.info,file="/home/memon/projects/msdrp/wikipath/hsaSPIA.RData")
#####################

load("/home/memon/projects/msdrp/keggx/hsaSPIA.RData")
load("/home/memon/projects/msdrp/data/real_kegg/hsaSPIA.RData")

#####________SPIA____________#####
##~~~Just for Testin purpose
library(SPIA)

load("/home/memon/projects/msdrp/data/disease.genes50.lfc.namedVec.HGNC.RData")

spia_wiki = list()
for (i in 1:length(lfc_hgnc)) {
  spia_wiki[[1]] = spia(de = lfc_hgnc[[1]], all = hgnc_all, data.dir="/home/memon/projects/msdrp/wikipath/",organism="hsa")
}

names(spia_wiki) = names(lfc_hgnc)
save(spia_wiki,file = "/home/memon/projects/msdrp/data/spia/spia_wikipath_custom_results.RData")
# save(spia_kegg,file = "/home/memon/projects/msdrp/data/spia/spia_kegg_custom_results.RData")

# load("/home/memon/projects/msdrp/data/disease.genes50.lfc.namedVec.RData")
spia_wiki_ad1 = spia(de = lfc_hgnc[[3]], all = hgnc_all, data.dir="/home/memon/projects/msdrp/forspia/test/",organism="hsa")
# spia_wiki_ad2 = spia(de = lfc_entrez[[3]], all = entrez_all, data.dir="/home/memon/projects/msdrp/data/real_kegg/",organism="hsa")
# spia_wiki_ad3 = spia(de = lfc_entrezID[[3]], all = entrezID_all, data.dir="/home/memon/projects/drug_target_prioritization/GCMap/data/spia_react/",organism="hsa")
