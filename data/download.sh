#!/usr/bin/env bash

######################
# MSDRP Downloader ###
#########################################################
# This script takes care of downloading all resources ###
# you'll need to run the MSDRP pipeline               ###
#########################################################

# Download the GRASP data, if it doesn't already exist
if [ ! -f GraspFullDataset2.zip ]; then
    wget https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip
else
	echo "Already downloaded GRASP"
fi

# Download the GWAS Catalog data, if it doesn't already exist
if [ ! -f gwas_catalog_v1.0-associations_e93_r2019-01-31.tsv ]; then
    wget -o gwas_catalog_v1.0-associations_e93_r2019-01-31.tsv https://www.ebi.ac.uk/gwas/api/search/downloads/full
else
    echo "Already downloaded GWAS Catalog"
fi

# Download the PHEWAS Catalog data, if it doesn't already exist
if [ ! -f phewas-catalog.csv.zip ]; then
    wget https://phewascatalog.org/files/phewas-catalog.csv.zip
else
    echo "Already downloaded PHEWAS Catalog"
fi

# Download the GWAS DB data, if it doesn't already exist
if [ ! -f gwasdb_20150819_snp_trait.gz ]; then
    wget ftp://147.8.193.36/GWASdb/gwasdb_20150819_snp_trait.gz
else
    echo "Already downloaded GWAS DB"
fi
