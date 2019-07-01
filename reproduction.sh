#!/usr/bin/env bash

# Halt the script on any error. See: https://stackoverflow.com/questions/19622198/what-does-set-e-mean-in-a-bash-script
set -e

####################
# 1. Preprocessing #
####################

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#preprocess-gwas-with-stopgap2
echo $"\n1.1 Running STOPGAP\n"
Rscript R/preprocessing/STOPGAP2_run.R data/ R/preprocessing/STOPGAP2_functions.R

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#prepare-mesh-to-efo-mappings
echo $"\n1.2 Mapping MeSH to EFO\n"
Rscript R/preprocessing/MeSH2EFO.R data/

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#download-open-targets-data
echo $"\n1.3 Retrieving DEGs\n"
Rscript R/preprocessing/RetrieveDEGs.R data/

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#download-and-preprocess-lincs
echo $"\n1.4 Retreiving drug response data\n"
Rscript R/preprocessing/RetrieveDrugResponseData.R data/

###############
# 2. Analysis #
###############

echo $"\n2.1 Calculating significant overlaps\n"
Rscript R/analysis/SignificantOverlap.R data/

echo $"\n2.2 Calculating Disease SPIAs\n"
Rscript R/analysis/DiseaseSPIA.R data/

echo $"\n2.3 Calculating Drug SPIAs\n"
Rscript R/analysis/DrugSPIA.R data/

echo $"\n2.4 Checking Distributions\n"
Rscript R/analysis/CheckDistribution.R data/

echo $"\n2.5 Checking Drug Combinations\n"
Rscript R/analysis/DrugCombination.R data/
