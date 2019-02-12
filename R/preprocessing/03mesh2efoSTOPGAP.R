#' 3rd script
#' summary:
#' Mapping of MeSH terms to EFO IDs 

#' 

library(data.table) #v 1.11.4
#source("https://bioconductor.org/biocLite.R")
#biocLite("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

setwd("/home/memon/projects/msdrp/")


#xm <- fread("gunzip -c ./ctd_gene_dis.txt.gz", header = TRUE, skip = 1)

# options
# permissive
max.gene.rank.min <- Inf
min.gene.score <- -Inf
sources.to.exclude <- c("OMIM", "Orphanet")
## strict
#max.gene.rank.min <- 3
#min.gene.score <- 1
#sources.to.exclude <- c("OMIM", "Orphanet")

# read data in
#stopgap <- fread("./data/stopgap.gene.mesh.txt")
load("./data/stopgap.gene.mesh.RData")
stopgap <- stopgap.gene.mesh
rm(stopgap.gene.mesh)

length(unique(stopgap$msh))
length(unique(stopgap$snp.gwas))
length(unique(stopgap$snp.ld))

# keep relevant columns
stopgap = data.table(stopgap)
stopgap <- unique(stopgap[, .(snp.ld,gene.v19, msh, pvalue, gene.score, gene.rank.min, source)])
# replace p-values of zero with arbitrarily low p-value
stopgap[pvalue == 0, pvalue := 3e-324]

# map gene symbols to Ensembl
stopgap.genes <- unique(stopgap[, gene.v19])
anno <- as.data.table(select(EnsDb.Hsapiens.v86, keys = stopgap.genes, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL")))
stopgap <- merge(stopgap, anno, by.x = "gene.v19", by.y = "SYMBOL", all = FALSE)

# map MeSH terms to EFO using Zooma 
mesh.terms <- gsub(" ", "+", unique(stopgap[, msh]))
#fwrite(data.frame(mesh.terms),"mesh.terms.txt",quote = F)

#write.table(mesh.terms,file = "/home/memon/projects/STOPGAP-master/STOPGAP_data/meshterm.txt",quote = F,col.names = F,row.names = F)
#. this api connection sometimes doesn't work # doesn't event work on my windows machine
zooma <- foreach(i = seq(mesh.terms), .combine = rbind) %dopar% {
    tmp <- fromJSON(content(GET(paste0("http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue=", mesh.terms[i]))))
    if (length(tmp) > 0) {
        tmp <- data.table(mesh.term = gsub("\\+", " ", mesh.terms[i]), efo.term = tmp$annotatedProperty$propertyValue[1], efo.url = tmp$semanticTags[[1]][1], confidence = tmp$confidence[1])
    }
}

# save(zooma,file = "./data/zooma_results.RData")
load("./data/zooma_results.RData")
# merge
stopgap <- merge(stopgap, zooma, by.x = "msh", by.y = "mesh.term", all = FALSE)
# filter
stopgap <- stopgap[gene.rank.min <= max.gene.rank.min & gene.score > min.gene.score & !(source %in% sources.to.exclude), ]
# tidy
stopgap[, efo.id := sub(".+\\/([A-Za-z]+_[0-9]+)", "\\1", efo.url)]
stopgap <- stopgap[, .(snp.ld,ensembl.id = GENEID, gene.symbol = gene.v19, efo.id, efo.term, stopgap.pvalue = pvalue, stopgap.gene.score = gene.score, stopgap.gene.rank = gene.rank.min)]
unique(stopgap$efo.id)
fwrite(stopgap, "./data/stopgap_snp.tsv", sep = "\t")


#zoom <- read.csv("/home/memon/projects/drug_target_prioritization/GCMap/data/zooma.tsv",header = T,sep = "\t")

length(unique(stopgap2$ensembl.id))
length(unique(stopgap2$gene.symbol))
length(unique(stopgap2$efo.id))
length(unique(stopgap2$snp.ld))
