#' 1st script
#' summary:
#' cleaning up STOPGAP data


library(data.table)
library(ggplot2)
library(riverplot)
library(RColorBrewer)
library(httr)
library(jsonlite)
library(foreach)
library(dplyr)
library(tidyr)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

setwd("/home/memon/projects/msdrp/")

# read efo.ids to therapeutic area mapping file
load("./data/drug2disease.therapeutic.RData")
therap.area <- unique(drug2disease.therapeutic[,c(1,3)])
length(unique(therap.area$efo.id))
length(unique(therap.area$therapeutic.area))
therap.area$therapeutic.area = trimws(therap.area$therapeutic.area)
opentargets.tas <- data.table(unique(therap.area))

relevant.tas <- c("cardiovascular disease", "digestive system disease", "endocrine system disease", "hematological system disease", "immune system disease", "infectious disease", "metabolic disease", "neoplasm", "nervous system disease", "reproductive system disease", "respiratory system disease", "skeletal system disease")

fisher.genes <- fread("./data/fisher.genes.tsv")
fisher.drugs <- fread("./data/fisher.drugs.tsv")

drug_dis = unique(fisher.drugs[existing.indication == TRUE & p.adjusted < 0.05])
drug_dis_new = unique(fisher.drugs[existing.indication == FALSE & p.adjusted < 1e-10])

length(unique(fisher.genes$efo.id.stopgap))


# efo.ids <- unique(fisher.genes$efo.id.stopgap)
## hypothesis 1 (fisher.genes)
# add therapeutic area
# fisher.genes <- merge(fisher.genes, opentargets.tas[, .(ID, therapeutic.area = TA_LABEL)], by.x = "efo.id.stopgap", by.y = "ID", all = FALSE)
fisher.genes <- merge(fisher.genes, opentargets.tas[, .(efo.id, therapeutic.area)], by.x = "efo.id.stopgap", by.y = "efo.id", all = FALSE)
fisher.genes <- unique(fisher.genes[,-4]) # give false counts if this column remains
fisher.genes <- unique(fisher.genes)
# calculate number of disease per TA
genes.numbers <- fisher.genes[same.disease == TRUE & therapeutic.area %in% relevant.tas, .(number = sum(p.adjusted < 0.05)), by = therapeutic.area]
# plot
png("./data/fisher.genes.ta.barplot.png", res = 150, width = 9 * 150, height = 6 * 150)
print(ggplot(genes.numbers, aes(x = reorder(therapeutic.area, -round(number)), y = round(number))) +
      geom_bar(stat = "identity", colour = "black", fill = "#0066ff") +
      xlab("Therapeutic area") +
      ylab("Number of significant diseases") +
      theme_bw(18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))
dev.off()

## hypothesis 2 (fisher.drugs)
# add therapeutic area
#fisher.drugs <- merge(fisher.drugs, opentargets.tas[, .(ID, therapeutic.area = TA_LABEL)], by.x = "efo.id", by.y = "ID", all = FALSE, allow.cartesian = TRUE)
fisher.drugs <- merge(fisher.drugs, opentargets.tas[, .(efo.id, therapeutic.area)], by.x = "efo.id", by.y = "efo.id", all = FALSE, allow.cartesian = TRUE)
# calculate number of drugs per TA
drugs.numbers <- unique(fisher.drugs[existing.indication == TRUE & therapeutic.area %in% relevant.tas, .(number = sum(p.adjusted < 0.05)), by = therapeutic.area])
# plot
png("./data/fisher.drugs.ta.barplot.png", res = 150, width = 9 * 150, height = 6 * 150)
print(ggplot(drugs.numbers, aes(x = reorder(therapeutic.area, -number), y = number)) +
      geom_bar(stat = "identity", colour = "black", fill = "#ff6600") +
      xlab("Therapeutic area") +
      ylab("Number of significant diseases") +
      theme_bw(18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))
dev.off()

## repurposing visualitation with sankey plot
# grab sources (repositioning from)

load("./data/drugGWAS.genes50.padj1e-5.RData")
length(unique(drugGWAS.genes$efo.id))
length(unique(opentargets.tas$efo.id))
# add therapeutic area
drugGWAS.genes <- merge(drugGWAS.genes, opentargets.tas[, .(efo.id, therapeutic.area)], by= "efo.id", all = FALSE, allow.cartesian = TRUE)

sources <- unique(drugGWAS.genes[therapeutic.area %in% relevant.tas & existing.indication == TRUE, .(chembl.id, from.ta = therapeutic.area)])
# grab targets (repurposing to)
targets <- unique(drugGWAS.genes[therapeutic.area %in% relevant.tas & existing.indication == FALSE & p.adjusted < 1e-10, .(chembl.id, to.ta = therapeutic.area)])
# merge
drugs.repos <- merge(sources, targets, by = "chembl.id", all = FALSE, allow.cartesian = TRUE)
# edges
edges <- drugs.repos[, .(Value = .N^2), by = .(N1 = paste(from.ta, 1), N2 = paste(to.ta, 2))]
# nodes
nodes <- data.table(ID = c(paste(relevant.tas, 1), paste(relevant.tas, 2)), x = c(rep(1, length(relevant.tas)), rep(2, length(relevant.tas))), labels = rep(relevant.tas, 2), col = rep(brewer.pal(length(relevant.tas), "Set3"), 2))
# create riverplot object
sankey <- makeRiver(as.data.frame(nodes), as.data.frame(edges))
png("./data/drugGWAS.genes.repositioning.sankeyplot.png", res = 600, width = 12 * 600, height = 30 * 600)
plot(sankey, srt = "0", textcex = 1.4)
dev.off()

## explore trajectories
# most common
edges[head(order(-Value), 10)]
# largest sources
edges[, .(Value = sum(Value)), by = N1][order(-Value),]
# largest targets
edges[, .(Value = sum(Value)), by = N2][order(-Value),]



