#' 1st script
#' summary:
#' Mapping of MeSH terms to EFO IDs from GWASs Data processed by STOPGAP Pipeline
#' for downloading Differential Gene Expression Data from Open Targtets by EFO IDs
#' EFO ontology (version 2.105) file has been used for these mapping purpose

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(biomaRt)))

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~: get GWASs Data from STOPGAP pipeline output :~~~~~~~~~~~~~~~~~~#

load(file.path(dataFolder, "stopgap.gene.mesh.RData"))
GWASs = data.table(stopgap.gene.mesh)
rm(stopgap.gene.mesh)

sprintf("Number of unique MeSH terms: %d", length(unique(GWASs$msh)))
sprintf("Number of unique genes: %d", length(unique(GWASs$gene.v19)))
sprintf("Number of unique SNPs: %d", length(unique(GWASs$snp.ld)))

# keep relevant columns and remove redundant rows
GWASs = unique(GWASs[, .(snp.ld, gene, msh, pvalue, gene.score, gene.rank.min, source)])

sprintf("Number of p-values with underflow: %d", nrow(GWASs[pvalue == 0]))
# replace p-values of zero with arbitrarily low p-value
GWASs[pvalue == 0, pvalue := 3e-324]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~: get mapping among ENTREZ, HGNC & ENSEMBL_IDs using Biomart Package :~~~~~~#

ensembl = useEnsembl(biomart="ensembl", version=97, dataset="hsapiens_gene_ensembl")
val = c(1:23,"X","Y")
gene_id = getBM(attributes=c('entrezgene_id','hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position'),
                filters ='chromosome_name', values =val, mart = ensembl)

rm(ensembl,val)

gene_id$entrezgene = gsub("^$", NA, gene_id$ensembl_gene_id)
gene_id = gene_id[which(! is.na(gene_id$ensembl_gene_id)),]
gene_id = data.table(gene_id)
gene_id = unique(gene_id[, c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')])
names(gene_id) = c("ENTREZ", "ensembl.id", "HGNC")
gene_id = gene_id[! duplicated(gene_id$ensembl.id),]

save(gene_id, file = file.path(dataFolder, "geneID_v97.RData"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~: map Gene Symbols in GWAS Data to Ensembl Identifiers :~~~~~~~~~~~~~~~#

sprintf("Number of unique GWAS genes before mapping: %d", length(unique(GWASs$gene)))

GWASs = merge(GWASs, gene_id, by.x = "gene", by.y = "HGNC", all = FALSE)
sprintf("Number of unique GWAS genes after mapping: %d", length(unique(GWASs$gene)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~: map MeSH terms in GWAS dataset to EFO IDs :~~~~~~~~~~~~~~~~~~~#

#' read EFO ontolgy (version 2.105) files in csv format downloaded from bioportal on 11th March, 2019
#' "http://data.bioontology.org/ontologies/EFO/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv"

efo = fread(file.path(dataFolder, "EFO.csv"))
efo = efo[, c(1, 2, 3)] # filter out unnecessary columns
efo$`Class ID` = gsub(".*efo/", "", efo$`Class ID`) # replace urls to get only EFO IDs
efo = efo[efo$`Class ID` %like% "EFO",] # filter out other disease IDs other than EFO
efo = unite(efo, Disease, c('Preferred Label', Synonyms), sep = "|", remove = FALSE)
efo = efo[, c(1, 2)]

# create new rows for each synonyms of EFO IDs 
efo = unique(efo %>%
    mutate(Disease = strsplit(as.character(Disease), "\\|")) %>%
    unnest(Disease))
efo$Disease = gsub("^$", NA, efo$Disease)
efo = efo[which(! is.na(efo$Disease)),]
efo$upper = toupper(efo$Disease)
efo$upper = trimws(efo$upper, which = "both")

mesh.terms = unique(GWASs[, 3])
mesh.terms$upper = toupper(mesh.terms$msh)
mesh.terms$upper = trimws(mesh.terms$upper, which = "both")

# merge MeSH terms from GWASs to efo IDs
mesh2efo = merge(mesh.terms, efo, by = "upper")
mesh2efo = mesh2efo[, c(2, 3, 4)]
names(mesh2efo) = c("mesh", "efo.id", "efo.term")
mesh2efo = mesh2efo[!duplicated(mesh2efo[,c('efo.id')]),]

GWASs = merge(GWASs, mesh2efo, by.x = "msh", by.y = "mesh", all = FALSE)
GWASs = GWASs[, .(snp.ld, ensembl.id = ensembl.id, gene.symbol = gene, efo.id, efo.term, GWASs.pvalue = pvalue, GWASs.gene.score = gene.score, GWASs.gene.rank = gene.rank.min, source)]
sprintf("Number of unique EFO IDs after final mapping: %d", length(unique(GWASs$efo.id)))
sprintf("Number of unique Ensembl IDs after final mapping: %d", length(unique(GWASs$ensembl.id)))

save(GWASs, file = file.path(dataFolder, "GWASs.RData"))
